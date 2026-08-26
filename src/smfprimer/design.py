"""Minimal, dependency-free primer candidate design."""

from __future__ import annotations

from .chemistry import (
    encode_converted_template,
    gc_fraction,
    normalize_sequence,
    reverse_complement,
    wallace_tm,
)
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


def _candidate(
    unconverted_template: str,
    encoded_template: str,
    start: int,
    end: int,
    side: str,
    parameters: DesignParameters,
) -> CandidatePrimer | None:
    unconverted_site = unconverted_template[start:end]
    if set(unconverted_site) - set("ACGT"):
        return None
    binding_site = encoded_template[start:end]
    sequence = binding_site if side == "forward" else reverse_complement(binding_site)
    unconverted_sequence = (
        unconverted_site if side == "forward" else reverse_complement(unconverted_site)
    )
    tm = wallace_tm(sequence)
    gc = gc_fraction(sequence)
    if not parameters.min_tm <= tm <= parameters.max_tm:
        return None
    if not parameters.min_gc <= gc <= parameters.max_gc:
        return None
    degeneracies = sequence.count("R") + sequence.count("Y")
    score = (
        degeneracies * 100.0
        + abs(len(sequence) - parameters.optimum_length) * 2.0
        + abs(tm - parameters.optimum_tm)
        + abs(gc - 0.5) * 10.0
    )
    return CandidatePrimer(
        sequence=sequence,
        unconverted_sequence=unconverted_sequence,
        side=side,
        start=start,
        end=end,
        tm=tm,
        gc_fraction=gc,
        degeneracies=degeneracies,
        score=score,
    )


def _candidates(
    unconverted_template: str,
    encoded_template: str,
    target_start: int,
    target_end: int,
    side: str,
    parameters: DesignParameters,
) -> list[CandidatePrimer]:
    if side == "forward":
        starts = range(max(0, target_start - parameters.search_window), target_start)
        sites = (
            (start, start + length)
            for start in starts
            for length in range(parameters.min_length, parameters.max_length + 1)
            if start + length <= target_start
        )
    else:
        ends = range(
            target_end + 1,
            min(len(encoded_template), target_end + parameters.search_window) + 1,
        )
        sites = (
            (end - length, end)
            for end in ends
            for length in range(parameters.min_length, parameters.max_length + 1)
            if end - length >= target_end
        )
    candidates = [
        candidate
        for start, end in sites
        if (
            candidate := _candidate(
                unconverted_template, encoded_template, start, end, side, parameters
            )
        )
        is not None
    ]
    return sorted(candidates, key=lambda item: (item.score, item.start, item.end))


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
    """Design ranked primer pairs around a zero-based, half-open target interval.

    Coordinates in results always refer to the supplied top-strand reference,
    even when the bottom strand is the converted strand.
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
    forwards = _candidates(oriented, encoded, oriented_start, oriented_end, "forward", parameters)
    reverses = _candidates(oriented, encoded, oriented_start, oriented_end, "reverse", parameters)

    pairs: list[PrimerPair] = []
    for forward in forwards:
        compatible_for_forward = 0
        for reverse in reverses:
            amplicon_length = reverse.end - forward.start
            if (
                parameters.min_amplicon_size is not None
                and amplicon_length < parameters.min_amplicon_size
            ):
                continue
            if (
                parameters.max_amplicon_size is not None
                and amplicon_length > parameters.max_amplicon_size
            ):
                continue
            compatible_for_forward += 1
            if converted_strand is ConvertedStrand.TOP:
                amplicon_start, amplicon_end = forward.start, reverse.end
                reported_forward, reported_reverse = forward, reverse
            else:
                amplicon_start = len(reference) - reverse.end
                amplicon_end = len(reference) - forward.start
                reported_forward = _remap(forward, len(reference))
                reported_reverse = _remap(reverse, len(reference))
            pairs.append(
                PrimerPair(
                    reported_forward,
                    reported_reverse,
                    amplicon_start,
                    amplicon_end,
                    converted_strand,
                    forward.score + reverse.score,
                )
            )
            # For a fixed forward primer, reverses are score-sorted. No pair
            # below the first max_results compatible reverses can enter the
            # global top max_results.
            if compatible_for_forward >= parameters.max_results:
                break
    return sorted(pairs, key=lambda pair: pair.score)[: parameters.max_results]


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
        message = ""
        if not pairs:
            message = "no primer pairs passed the length, Tm, GC, and amplicon constraints"
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
