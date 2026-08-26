"""Evaluation of user-supplied primer pairs."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, replace
from pathlib import Path

from .chemistry import (
    encode_converted_template,
    normalize_sequence,
    reverse_complement,
)
from .models import (
    CandidatePrimer,
    ConvertedStrand,
    DesignOutcome,
    DesignParameters,
    DesignTarget,
    PrimerPair,
    PrimerSpecificity,
    TargetContext,
    Workflow,
)
from .primer3_engine import check_pair, problems
from .specificity import assess_specificity

_BASES = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("CG"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}


@dataclass(frozen=True)
class PrimerSet:
    """One named forward/reverse ordered-primer pair."""

    name: str
    forward: str
    reverse: str

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("primer-set name must not be empty")


@dataclass(frozen=True)
class PrimerMetrics:
    sequence: str
    length: int
    tm: float
    gc_fraction: float
    degeneracies: int
    score: float
    warnings: tuple[str, ...]

    @property
    def passes_constraints(self) -> bool:
        return not self.warnings


@dataclass(frozen=True)
class TemplateAmplicon:
    reference_name: str
    start: int
    end: int
    forward_start: int
    forward_end: int
    reverse_start: int
    reverse_end: int
    within_product_size: bool

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class PrimerPairEvaluation:
    name: str
    forward: PrimerMetrics
    reverse: PrimerMetrics
    pair_score: float
    template_amplicons: tuple[TemplateAmplicon, ...]
    warnings: tuple[str, ...]
    specificity: PrimerSpecificity | None = None

    @property
    def valid_template_amplicons(self) -> tuple[TemplateAmplicon, ...]:
        return tuple(product for product in self.template_amplicons if product.within_product_size)

    @property
    def status(self) -> str:
        if self.forward.warnings or self.reverse.warnings:
            return "constraint_warning"
        if not self.valid_template_amplicons:
            return "no_template_product"
        if self.specificity is not None and self.specificity.off_target_amplicon_count:
            return "off_targets"
        return "pass"


def evaluate_primer_pairs(
    primer_sets: PrimerSet | Iterable[PrimerSet],
    template: str,
    *,
    template_name: str = "template",
    workflow: Workflow | str = Workflow.CONVERSION,
    converted_strand: ConvertedStrand | str = ConvertedStrand.TOP,
    context: TargetContext | str = TargetContext.BOTH,
    parameters: DesignParameters | None = None,
    bowtie_index: str | Path | None = None,
    mismatches: int = 2,
    bowtie: str = "bowtie",
    maximum_expansions: int = 256,
) -> list[PrimerPairEvaluation]:
    """Score one or more ordered primer pairs against a reference template.

    The template is the unconverted top-strand sequence. IUPAC-compatible
    binding sites are located on its predicted converted strand. Supplying a
    Bowtie 1 index adds individual-hit and off-target-product metrics.
    """
    sets = [primer_sets] if isinstance(primer_sets, PrimerSet) else list(primer_sets)
    if any(not isinstance(primer_set, PrimerSet) for primer_set in sets):
        raise TypeError("primer_sets must contain PrimerSet instances")
    if not template_name:
        raise ValueError("template_name must not be empty")

    reference = normalize_sequence(template)
    workflow = Workflow(workflow)
    converted_strand = ConvertedStrand(converted_strand)
    context = TargetContext(context)
    parameters = parameters or DesignParameters()
    if converted_strand is ConvertedStrand.TOP:
        oriented = reference
    else:
        oriented = reverse_complement(reference)
    encoded = encode_converted_template(oriented, workflow, context)

    evaluations = [
        _evaluate_one(
            primer_set,
            reference,
            encoded,
            template_name,
            converted_strand,
            parameters,
        )
        for primer_set in sets
    ]
    if bowtie_index is None or not evaluations:
        return evaluations

    adapters = tuple(
        _specificity_adapter(evaluation, converted_strand) for evaluation in evaluations
    )
    target = DesignTarget.create(
        target_id=template_name,
        reference_name=template_name,
        sequence=reference,
        target_start=0,
        target_end=1,
        source="primer_evaluation",
    )
    outcome = DesignOutcome(
        target=target,
        workflow=workflow,
        context=context,
        converted_strand=converted_strand,
        parameters=parameters,
        pairs=adapters,
    )
    assessed = assess_specificity(
        [outcome],
        bowtie_index,
        mismatches=mismatches,
        bowtie=bowtie,
        maximum_expansions=maximum_expansions,
    )[0]
    return [
        replace(evaluation, specificity=pair.specificity)
        for evaluation, pair in zip(evaluations, assessed.pairs, strict=True)
    ]


def _evaluate_one(
    primer_set: PrimerSet,
    reference: str,
    encoded: str,
    template_name: str,
    converted_strand: ConvertedStrand,
    parameters: DesignParameters,
) -> PrimerPairEvaluation:
    primer3_result = check_pair(primer_set.forward, primer_set.reverse, parameters)
    forward = _metrics(primer_set.forward, primer3_result, "LEFT")
    reverse = _metrics(primer_set.reverse, primer3_result, "RIGHT")
    forward_sites = _binding_sites(forward.sequence, encoded)
    reverse_sites = _binding_sites(reverse_complement(reverse.sequence), encoded)
    products = _template_products(
        forward_sites,
        reverse_sites,
        len(reference),
        template_name,
        converted_strand,
        parameters,
    )

    warnings = []
    if not forward_sites:
        warnings.append("forward primer has no compatible template binding site")
    if not reverse_sites:
        warnings.append("reverse primer has no compatible template binding site")
    if forward_sites and reverse_sites and not products:
        warnings.append("primer sites do not form an inward-facing template product")
    elif products and not any(product.within_product_size for product in products):
        warnings.append("template products fall outside the configured product-size range")
    if sum(product.within_product_size for product in products) > 1:
        warnings.append("multiple template products satisfy the product-size range")
    return PrimerPairEvaluation(
        name=primer_set.name,
        forward=forward,
        reverse=reverse,
        pair_score=float(primer3_result["PRIMER_PAIR_0_PENALTY"]),
        template_amplicons=products,
        warnings=tuple(warnings),
    )


def _metrics(sequence: str, result: dict, side: str) -> PrimerMetrics:
    sequence = normalize_sequence(sequence)
    length = len(sequence)
    tm = float(result[f"PRIMER_{side}_0_TM"])
    gc = float(result[f"PRIMER_{side}_0_GC_PERCENT"]) / 100.0
    degeneracies = sum(base not in "ACGT" for base in sequence)
    return PrimerMetrics(
        sequence=sequence,
        length=length,
        tm=tm,
        gc_fraction=gc,
        degeneracies=degeneracies,
        score=float(result[f"PRIMER_{side}_0_PENALTY"]),
        warnings=problems(result, side),
    )


def _binding_sites(pattern: str, template: str) -> tuple[tuple[int, int], ...]:
    length = len(pattern)
    return tuple(
        (start, start + length)
        for start in range(len(template) - length + 1)
        if all(
            _BASES[primer_base] & _BASES[template_base]
            for primer_base, template_base in zip(
                pattern, template[start : start + length], strict=True
            )
        )
    )


def _template_products(
    forward_sites: tuple[tuple[int, int], ...],
    reverse_sites: tuple[tuple[int, int], ...],
    reference_length: int,
    template_name: str,
    converted_strand: ConvertedStrand,
    parameters: DesignParameters,
) -> tuple[TemplateAmplicon, ...]:
    products = []
    for forward_start, forward_end in forward_sites:
        for reverse_start, reverse_end in reverse_sites:
            if reverse_start <= forward_start:
                continue
            reported_forward_start, reported_forward_end = forward_start, forward_end
            reported_reverse_start, reported_reverse_end = reverse_start, reverse_end
            start, end = forward_start, reverse_end
            if converted_strand is ConvertedStrand.BOTTOM:
                start, end = reference_length - end, reference_length - start
                reported_forward_start, reported_forward_end = (
                    reference_length - forward_end,
                    reference_length - forward_start,
                )
                reported_reverse_start, reported_reverse_end = (
                    reference_length - reverse_end,
                    reference_length - reverse_start,
                )
            length = end - start
            within = (
                parameters.min_amplicon_size is None or length >= parameters.min_amplicon_size
            ) and (parameters.max_amplicon_size is None or length <= parameters.max_amplicon_size)
            products.append(
                TemplateAmplicon(
                    reference_name=template_name,
                    start=start,
                    end=end,
                    forward_start=reported_forward_start,
                    forward_end=reported_forward_end,
                    reverse_start=reported_reverse_start,
                    reverse_end=reported_reverse_end,
                    within_product_size=within,
                )
            )
    return tuple(sorted(products, key=lambda product: (product.start, product.end)))


def _specificity_adapter(
    evaluation: PrimerPairEvaluation,
    converted_strand: ConvertedStrand,
) -> PrimerPair:
    valid = evaluation.valid_template_amplicons
    product = valid[0] if valid else None
    forward_start = product.forward_start if product else 0
    reverse_start = product.reverse_start if product else 0
    forward = CandidatePrimer(
        sequence=evaluation.forward.sequence,
        unconverted_sequence=evaluation.forward.sequence,
        side="forward",
        start=forward_start,
        end=forward_start + evaluation.forward.length,
        tm=evaluation.forward.tm,
        gc_fraction=evaluation.forward.gc_fraction,
        degeneracies=evaluation.forward.degeneracies,
        score=evaluation.forward.score,
    )
    reverse = CandidatePrimer(
        sequence=evaluation.reverse.sequence,
        unconverted_sequence=evaluation.reverse.sequence,
        side="reverse",
        start=reverse_start,
        end=reverse_start + evaluation.reverse.length,
        tm=evaluation.reverse.tm,
        gc_fraction=evaluation.reverse.gc_fraction,
        degeneracies=evaluation.reverse.degeneracies,
        score=evaluation.reverse.score,
    )
    return PrimerPair(
        forward=forward,
        reverse=reverse,
        amplicon_start=product.start if product else 0,
        amplicon_end=product.end if product else 1,
        converted_strand=converted_strand,
        score=evaluation.pair_score,
    )
