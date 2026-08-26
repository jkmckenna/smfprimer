"""GenBank serialization for annotated primer-design records."""

from __future__ import annotations

import json
import re
import textwrap
from datetime import date

from .models import ConvertedStrand, DesignOutcome, PrimerPair, ReferenceFeature


def format_genbank(outcomes: list[DesignOutcome]) -> str:
    """Serialize design templates and ranked primer pairs as GenBank records."""
    return "".join(_record(outcome, index) for index, outcome in enumerate(outcomes, 1))


def _record(outcome: DesignOutcome, index: int) -> str:
    target = outcome.target
    sequence = target.sequence.upper()
    locus = _locus_name(target.target_id, index)
    today = date.today().strftime("%d-%b-%Y").upper()
    lines = [
        f"LOCUS       {locus:<16} {len(sequence):>11} bp    DNA     "
        f"{target.topology:<8} UNA {today}",
        f"DEFINITION  smfprimer design template for {target.target_id}.",
        f"ACCESSION   {locus}",
        f"VERSION     {locus}",
        "KEYWORDS    primer design; single-molecule footprinting.",
        f"SOURCE      {target.reference_name}",
        "  ORGANISM  synthetic construct",
        "FEATURES             Location/Qualifiers",
    ]
    lines.extend(
        _feature(
            "source",
            _location(0, len(sequence)),
            [
                ("label", target.target_id),
                ("mol_type", "other DNA"),
                ("reference_name", target.reference_name),
                ("reference_start", target.sequence_start),
                ("source_mode", target.source),
            ],
        )
    )
    for feature in target.reference_features:
        lines.extend(_reference_feature(feature))
    lines.extend(
        _feature(
            "misc_feature",
            _location(target.target_start, target.target_end),
            [
                ("label", f"{target.target_id}_required_target"),
                ("note", "required target interval"),
                ("reference_start", target.reference_target_start),
                ("reference_end", target.reference_target_end),
            ],
        )
    )
    for rank, pair in enumerate(outcome.pairs, 1):
        lines.extend(_pair_features(outcome, pair, rank))
    lines.append("ORIGIN")
    lines.extend(_origin(sequence))
    lines.append("//")
    return "\n".join(lines) + "\n"


def _reference_feature(feature: ReferenceFeature) -> list[str]:
    qualifiers = [(name, value) for name, values in feature.qualifiers for value in values]
    if feature.label and not any(name == "label" for name, _ in qualifiers):
        qualifiers.insert(0, ("label", feature.label))
    return _feature(
        feature.type,
        _segments_location(feature),
        qualifiers,
    )


def _segments_location(feature: ReferenceFeature) -> str:
    parts = [_location(segment.start, segment.end) for segment in feature.segments]
    location = parts[0] if len(parts) == 1 else f"join({','.join(parts)})"
    return f"complement({location})" if feature.strand == -1 else location


def _pair_features(
    outcome: DesignOutcome,
    pair: PrimerPair,
    rank: int,
) -> list[str]:
    target = outcome.target
    specificity = pair.specificity
    amplicon_qualifiers: list[tuple[str, object]] = [
        ("label", f"rank_{rank}_amplicon"),
        ("rank", rank),
        ("amplicon_length", pair.amplicon_length),
        ("pair_score", round(pair.score, 3)),
        ("design_engine", "primer3"),
        ("workflow", outcome.workflow.value),
        ("converted_strand", outcome.converted_strand.value),
        ("reference_start", target.sequence_start + pair.amplicon_start),
        ("reference_end", target.sequence_start + pair.amplicon_end),
    ]
    if specificity is not None:
        amplicon_qualifiers.extend(
            [
                ("specificity_status", specificity.status),
                ("specificity_index", specificity.index),
                ("specificity_mismatches", specificity.mismatches),
                ("specificity_forward_hits", specificity.forward_hit_count),
                ("specificity_reverse_hits", specificity.reverse_hit_count),
                (
                    "specificity_off_target_amplicons",
                    specificity.off_target_amplicon_count,
                ),
            ]
        )

    lines = _feature(
        "misc_feature",
        _location(pair.amplicon_start, pair.amplicon_end),
        amplicon_qualifiers,
    )
    forward_reverse = outcome.converted_strand is ConvertedStrand.BOTTOM
    reverse_reverse = outcome.converted_strand is ConvertedStrand.TOP
    lines.extend(
        _primer_feature(
            outcome,
            pair,
            rank,
            "forward",
            reverse=forward_reverse,
        )
    )
    lines.extend(
        _primer_feature(
            outcome,
            pair,
            rank,
            "reverse",
            reverse=reverse_reverse,
        )
    )
    return lines


def _primer_feature(
    outcome: DesignOutcome,
    pair: PrimerPair,
    rank: int,
    side: str,
    *,
    reverse: bool,
) -> list[str]:
    primer = pair.forward if side == "forward" else pair.reverse
    target = outcome.target
    location = _location(primer.start, primer.end, reverse=reverse)
    return _feature(
        "primer_bind",
        location,
        [
            ("label", f"rank_{rank}_{side}_primer"),
            ("rank", rank),
            ("primer_side", side),
            ("order_sequence", primer.sequence),
            ("unconverted_sequence", primer.unconverted_sequence),
            ("estimated_tm", round(primer.tm, 3)),
            ("thermodynamic_model", "primer3"),
            ("gc_fraction", round(primer.gc_fraction, 5)),
            ("degeneracies", primer.degeneracies),
            ("primer_score", round(primer.score, 3)),
            ("reference_start", target.sequence_start + primer.start),
            ("reference_end", target.sequence_start + primer.end),
        ],
    )


def _feature(
    key: str,
    location: str,
    qualifiers: list[tuple[str, object]],
) -> list[str]:
    lines = [f"     {key[:15]:<16}{location}"]
    for name, value in qualifiers:
        lines.extend(_qualifier(name, value))
    return lines


def _qualifier(name: str, value: object) -> list[str]:
    if isinstance(value, (dict, list, tuple)):
        text = json.dumps(value, sort_keys=True)
    else:
        text = str(value)
    text = text.replace('"', "'")
    prefix = f'                     /{name}="'
    if not any(character.isspace() for character in text):
        return [prefix + text + '"']
    wrapped = textwrap.wrap(
        text,
        width=max(1, 79 - len(prefix)),
        break_long_words=False,
        break_on_hyphens=False,
    ) or [""]
    if len(wrapped) == 1:
        return [prefix + wrapped[0] + '"']
    lines = [prefix + wrapped[0]]
    lines.extend("                     " + part for part in wrapped[1:-1])
    lines.append("                     " + wrapped[-1] + '"')
    return lines


def _location(start: int, end: int, *, reverse: bool = False) -> str:
    if not 0 <= start < end:
        raise ValueError("GenBank feature must satisfy 0 <= start < end")
    location = str(start + 1) if end == start + 1 else f"{start + 1}..{end}"
    return f"complement({location})" if reverse else location


def _origin(sequence: str) -> list[str]:
    lines = []
    sequence = sequence.lower()
    for start in range(0, len(sequence), 60):
        chunk = sequence[start : start + 60]
        groups = " ".join(chunk[index : index + 10] for index in range(0, len(chunk), 10))
        lines.append(f"{start + 1:>9} {groups}")
    return lines


def _locus_name(target_id: str, index: int) -> str:
    sanitized = re.sub(r"[^A-Za-z0-9_.-]", "_", target_id).strip("_.-") or "target"
    suffix = f"_{index}"
    return sanitized[: 16 - len(suffix)] + suffix
