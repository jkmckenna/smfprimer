"""Conversion-aware primer design powered by Primer3."""

from __future__ import annotations

from .chemistry import encode_converted_template, normalize_sequence, reverse_complement
from .models import (
    CandidatePrimer,
    ConvertedStrand,
    DesignOutcome,
    DesignParameters,
    DesignTarget,
    PrimerPair,
    TargetContext,
    Workflow,
)
from .primer3_engine import design


def design_primers(
    reference: str,
    target_start: int,
    target_end: int,
    *,
    workflow: Workflow | str = Workflow.CONVERSION,
    converted_strand: ConvertedStrand | str = ConvertedStrand.TOP,
    context: TargetContext | str = TargetContext.BOTH,
    parameters: DesignParameters | None = None,
) -> list[PrimerPair]:
    """Design Primer3-ranked pairs around a zero-based, half-open target.

    Coordinates in results always refer to the supplied top-strand reference,
    even when the bottom strand is the converted strand. Primer3 chooses sites
    on the concrete converted allele; the reported oligos are reconstructed
    from smfprimer's IUPAC-encoded template so conversion ambiguity is retained.
    """
    reference = normalize_sequence(reference)
    workflow = Workflow(workflow)
    converted_strand = ConvertedStrand(converted_strand)
    context = TargetContext(context)
    parameters = parameters or DesignParameters()
    if not 0 <= target_start < target_end <= len(reference):
        raise ValueError("target must satisfy 0 <= start < end <= sequence length")

    if converted_strand is ConvertedStrand.TOP:
        oriented = reference
        oriented_start, oriented_end = target_start, target_end
    else:
        oriented = reverse_complement(reference)
        oriented_start, oriented_end = len(reference) - target_end, len(reference) - target_start

    encoded = encode_converted_template(oriented, workflow, context)
    window_start = max(0, oriented_start - parameters.search_window)
    window_end = min(len(oriented), oriented_end + parameters.search_window)
    window_encoded = encoded[window_start:window_end]
    window_unconverted = oriented[window_start:window_end]
    local_start = oriented_start - window_start
    local_end = oriented_end - window_start

    result = design(window_encoded, local_start, local_end, parameters)
    pairs = []
    for index in range(int(result.get("PRIMER_PAIR_NUM_RETURNED", 0))):
        left_start, left_length = result[f"PRIMER_LEFT_{index}"]
        right_three_prime, right_length = result[f"PRIMER_RIGHT_{index}"]
        right_start = right_three_prime - right_length + 1
        left = _primer_from_result(
            result,
            index,
            "LEFT",
            "forward",
            left_start,
            left_length,
            window_start,
            window_encoded,
            window_unconverted,
        )
        right = _primer_from_result(
            result,
            index,
            "RIGHT",
            "reverse",
            right_start,
            right_length,
            window_start,
            window_encoded,
            window_unconverted,
        )
        if parameters.max_degeneracies is not None and (
            left.degeneracies > parameters.max_degeneracies
            or right.degeneracies > parameters.max_degeneracies
        ):
            continue
        amplicon_start = left.start
        amplicon_end = right.end
        if converted_strand is ConvertedStrand.BOTTOM:
            amplicon_start = len(reference) - right.end
            amplicon_end = len(reference) - left.start
            left = _remap(left, len(reference))
            right = _remap(right, len(reference))
        pairs.append(
            PrimerPair(
                forward=left,
                reverse=right,
                amplicon_start=amplicon_start,
                amplicon_end=amplicon_end,
                converted_strand=converted_strand,
                score=float(result[f"PRIMER_PAIR_{index}_PENALTY"]),
            )
        )
    return pairs


def _primer_from_result(
    result: dict,
    index: int,
    primer3_side: str,
    side: str,
    local_start: int,
    length: int,
    window_start: int,
    encoded: str,
    unconverted: str,
) -> CandidatePrimer:
    local_end = local_start + length
    encoded_site = encoded[local_start:local_end]
    unconverted_site = unconverted[local_start:local_end]
    if side == "reverse":
        sequence = reverse_complement(encoded_site)
        unconverted_sequence = reverse_complement(unconverted_site)
    else:
        sequence = encoded_site
        unconverted_sequence = unconverted_site
    return CandidatePrimer(
        sequence=sequence,
        unconverted_sequence=unconverted_sequence,
        side=side,
        start=window_start + local_start,
        end=window_start + local_end,
        tm=float(result[f"PRIMER_{primer3_side}_{index}_TM"]),
        gc_fraction=float(result[f"PRIMER_{primer3_side}_{index}_GC_PERCENT"]) / 100.0,
        degeneracies=sum(base not in "ACGT" for base in sequence),
        score=float(result[f"PRIMER_{primer3_side}_{index}_PENALTY"]),
    )


def design_targets(
    targets: list[DesignTarget],
    *,
    workflow: Workflow | str = Workflow.CONVERSION,
    converted_strand: ConvertedStrand | str = ConvertedStrand.TOP,
    context: TargetContext | str = TargetContext.BOTH,
    parameters: DesignParameters | None = None,
) -> list[DesignOutcome]:
    """Design primers for normalized targets while retaining failed targets."""
    workflow = Workflow(workflow)
    converted_strand = ConvertedStrand(converted_strand)
    context = TargetContext(context)
    parameters = parameters or DesignParameters()
    outcomes = []
    for target in targets:
        pairs = design_primers(
            target.sequence,
            target.target_start,
            target.target_end,
            workflow=workflow,
            converted_strand=converted_strand,
            context=context,
            parameters=parameters,
        )
        message = "" if pairs else "Primer3 returned no primer pairs for the configured constraints"
        outcomes.append(
            DesignOutcome(
                target=target,
                workflow=workflow,
                context=context,
                converted_strand=converted_strand,
                parameters=parameters,
                pairs=tuple(pairs),
                message=message,
            )
        )
    return outcomes


def _remap(primer: CandidatePrimer, reference_length: int) -> CandidatePrimer:
    return CandidatePrimer(
        sequence=primer.sequence,
        unconverted_sequence=primer.unconverted_sequence,
        side=primer.side,
        start=reference_length - primer.end,
        end=reference_length - primer.start,
        tm=primer.tm,
        gc_fraction=primer.gc_fraction,
        degeneracies=primer.degeneracies,
        score=primer.score,
    )
