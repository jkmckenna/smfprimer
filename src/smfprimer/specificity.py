"""Primer-pair specificity assessment using Bowtie 1."""

from __future__ import annotations

import subprocess
import tempfile
from collections import defaultdict
from dataclasses import dataclass, replace
from pathlib import Path

from .models import DesignOutcome, OffTargetAmplicon, PrimerSpecificity

_IUPAC = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "S": "CG",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}


@dataclass(frozen=True)
class _Hit:
    reference_name: str
    start: int
    end: int
    strand: str
    mismatches: int


def expand_degenerate(sequence: str, *, maximum: int = 256) -> tuple[str, ...]:
    """Expand an IUPAC sequence into concrete DNA sequences."""
    expansions = [""]
    for base in sequence.upper():
        try:
            choices = _IUPAC[base]
        except KeyError as error:
            raise ValueError(
                f"unsupported primer base for specificity assessment: {base}"
            ) from error
        if len(expansions) * len(choices) > maximum:
            raise ValueError(
                f"primer {sequence} exceeds the specificity expansion limit of {maximum}"
            )
        expansions = [prefix + choice for prefix in expansions for choice in choices]
    return tuple(expansions)


def assess_specificity(
    outcomes: list[DesignOutcome],
    index: str | Path,
    *,
    mismatches: int = 2,
    bowtie: str = "bowtie",
    maximum_expansions: int = 256,
) -> list[DesignOutcome]:
    """Annotate designed pairs with products supported by Bowtie alignments.

    The Bowtie index defines the exact background being assessed. Any product
    that does not exactly match the designed reference coordinates is reported
    as an off-target product.
    """
    if not 0 <= mismatches <= 3:
        raise ValueError("Bowtie mismatch count must be between 0 and 3")
    if maximum_expansions < 1:
        raise ValueError("specificity expansion limit must be positive")

    query_map: dict[str, tuple[int, int, str]] = {}
    records: list[tuple[str, str]] = []
    for outcome_index, outcome in enumerate(outcomes):
        for pair_index, pair in enumerate(outcome.pairs):
            for side, primer in (("forward", pair.forward), ("reverse", pair.reverse)):
                for expansion_index, sequence in enumerate(
                    expand_degenerate(primer.sequence, maximum=maximum_expansions)
                ):
                    name = f"o{outcome_index}_p{pair_index}_{side}_e{expansion_index}"
                    query_map[name] = (outcome_index, pair_index, side)
                    records.append((name, sequence))

    if not records:
        return outcomes

    with tempfile.TemporaryDirectory(prefix="smfprimer-bowtie-") as directory:
        query_path = Path(directory) / "primers.fa"
        sam_path = Path(directory) / "alignments.sam"
        query_path.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in records))
        command = [
            bowtie,
            "-f",
            "-v",
            str(mismatches),
            "-a",
            "--best",
            "-S",
            "--sam-nohead",
            str(index),
            str(query_path),
            str(sam_path),
        ]
        try:
            completed = subprocess.run(command, capture_output=True, text=True, check=False)
        except OSError as error:
            raise ValueError(f"could not execute Bowtie 1 ({bowtie}): {error}") from error
        if completed.returncode:
            detail = completed.stderr.strip() or completed.stdout.strip() or "unknown error"
            raise ValueError(f"Bowtie 1 specificity assessment failed: {detail}")
        hits = _read_sam(sam_path, query_map)

    annotated: list[DesignOutcome] = []
    for outcome_index, outcome in enumerate(outcomes):
        pairs = []
        for pair_index, pair in enumerate(outcome.pairs):
            pair_hits = hits.get((outcome_index, pair_index), {})
            forward_hits = tuple(pair_hits.get("forward", {}).values())
            reverse_hits = tuple(pair_hits.get("reverse", {}).values())
            products = _products(forward_hits, reverse_hits, outcome.parameters)
            expected = (
                outcome.target.reference_name,
                outcome.target.sequence_start + pair.amplicon_start,
                outcome.target.sequence_start + pair.amplicon_end,
            )
            intended = sum(
                (product.reference_name, product.start, product.end) == expected
                for product in products
            )
            off_targets = tuple(
                product
                for product in products
                if (product.reference_name, product.start, product.end) != expected
            )
            specificity = PrimerSpecificity(
                index=str(index),
                mismatches=mismatches,
                forward_hit_count=len(forward_hits),
                reverse_hit_count=len(reverse_hits),
                intended_amplicon_count=intended,
                off_target_amplicons=off_targets,
            )
            pairs.append(replace(pair, specificity=specificity))
        annotated.append(replace(outcome, pairs=tuple(pairs)))
    return annotated


def _read_sam(
    path: Path,
    query_map: dict[str, tuple[int, int, str]],
) -> dict[tuple[int, int], dict[str, dict[tuple[str, int, int, str], _Hit]]]:
    hits: dict[tuple[int, int], dict[str, dict[tuple[str, int, int, str], _Hit]]] = defaultdict(
        lambda: defaultdict(dict)
    )
    with path.open() as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.strip() or line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                raise ValueError(f"invalid Bowtie SAM record at line {line_number}")
            name, flag_text, reference_name, position_text = fields[:4]
            if name not in query_map:
                raise ValueError(f"Bowtie returned an unknown query name: {name}")
            flag = int(flag_text)
            if flag & 4:
                continue
            start = int(position_text) - 1
            length = len(fields[9])
            strand = "-" if flag & 16 else "+"
            mismatches = next(
                (int(field[5:]) for field in fields[11:] if field.startswith("NM:i:")), 0
            )
            hit = _Hit(reference_name, start, start + length, strand, mismatches)
            outcome_index, pair_index, side = query_map[name]
            key = (reference_name, hit.start, hit.end, strand)
            current = hits[(outcome_index, pair_index)][side].get(key)
            if current is None or hit.mismatches < current.mismatches:
                hits[(outcome_index, pair_index)][side][key] = hit
    return hits


def _products(forward_hits, reverse_hits, parameters) -> tuple[OffTargetAmplicon, ...]:
    minimum = parameters.min_amplicon_size or 1
    maximum = parameters.max_amplicon_size
    products: dict[tuple[str, int, int], OffTargetAmplicon] = {}
    orientations = (
        (
            (hit for hit in forward_hits if hit.strand == "+"),
            (hit for hit in reverse_hits if hit.strand == "-"),
            "+",
        ),
        (
            (hit for hit in reverse_hits if hit.strand == "+"),
            (hit for hit in forward_hits if hit.strand == "-"),
            "-",
        ),
    )
    for left_hits, right_hits, forward_strand in orientations:
        right_by_reference: dict[str, list[_Hit]] = defaultdict(list)
        for hit in right_hits:
            right_by_reference[hit.reference_name].append(hit)
        for left in left_hits:
            for right in right_by_reference[left.reference_name]:
                length = right.end - left.start
                if right.start <= left.start or length < minimum:
                    continue
                if maximum is not None and length > maximum:
                    continue
                if forward_strand == "+":
                    forward_mismatches, reverse_mismatches = (
                        left.mismatches,
                        right.mismatches,
                    )
                else:
                    forward_mismatches, reverse_mismatches = (
                        right.mismatches,
                        left.mismatches,
                    )
                product = OffTargetAmplicon(
                    reference_name=left.reference_name,
                    start=left.start,
                    end=right.end,
                    forward_strand=forward_strand,
                    forward_mismatches=forward_mismatches,
                    reverse_mismatches=reverse_mismatches,
                )
                key = (product.reference_name, product.start, product.end)
                current = products.get(key)
                if current is None or (
                    product.forward_mismatches + product.reverse_mismatches
                    < current.forward_mismatches + current.reverse_mismatches
                ):
                    products[key] = product
    return tuple(
        sorted(products.values(), key=lambda item: (item.reference_name, item.start, item.end))
    )
