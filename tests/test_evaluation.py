from dataclasses import replace

import pytest

from smfprimer import (
    DesignParameters,
    PrimerSet,
    PrimerSpecificity,
    Workflow,
    evaluate_primer_pairs,
)
from smfprimer.chemistry import reverse_complement

FORWARD_SITE = "GAGAGATCTGGCAGCGGAGAG"
FORWARD_ORDER = "GAGAGATYTGGYAGYGGAGAG"
REVERSE_UNCONVERTED = "CTTTTCTGTCACCAATCCTGTCCC"
REVERSE_ORDER = "CTTTTCTRTCACCAATCCTRTCCC"
TEMPLATE = FORWARD_SITE + "A" * 20 + reverse_complement(REVERSE_UNCONVERTED)


def test_evaluate_one_known_deaminase_pair() -> None:
    result = evaluate_primer_pairs(
        PrimerSet("JM351_JM349", FORWARD_ORDER, REVERSE_ORDER),
        TEMPLATE,
        workflow=Workflow.DEAMINASE,
        parameters=DesignParameters(max_tm=72, min_amplicon_size=60, max_amplicon_size=70),
    )[0]
    assert result.name == "JM351_JM349"
    assert result.status == "pass"
    assert result.forward.length == 21
    assert result.forward.tm == pytest.approx(54.866, abs=0.001)
    assert result.forward.degeneracies == 3
    assert result.reverse.length == 24
    assert result.reverse.tm == pytest.approx(56.619, abs=0.001)
    assert result.reverse.degeneracies == 2
    assert result.pair_score == result.forward.score + result.reverse.score
    assert len(result.valid_template_amplicons) == 1
    assert result.valid_template_amplicons[0].start == 0
    assert result.valid_template_amplicons[0].end == len(TEMPLATE)


def test_evaluate_list_reports_constraint_and_product_warnings() -> None:
    results = evaluate_primer_pairs(
        [
            PrimerSet("known", FORWARD_ORDER, REVERSE_ORDER),
            PrimerSet("absent", "ACGTACGTACGTACGTAC", "TGCATGCATGCATGCATG"),
        ],
        TEMPLATE,
        workflow=Workflow.DEAMINASE,
        parameters=DesignParameters(
            min_tm=60,
            optimum_tm=62,
            min_amplicon_size=1,
            max_amplicon_size=50,
        ),
    )
    assert len(results) == 2
    assert results[0].status == "constraint_warning"
    assert "Temperature too low" in results[0].reverse.warnings
    assert "product-size" in results[0].warnings[0]
    assert results[1].status == "constraint_warning"
    assert any("binding site" in warning for warning in results[1].warnings)


def test_evaluate_bottom_strand_reports_top_reference_coordinates() -> None:
    results = evaluate_primer_pairs(
        PrimerSet("bottom", FORWARD_ORDER, REVERSE_ORDER),
        reverse_complement(TEMPLATE),
        workflow=Workflow.DEAMINASE,
        converted_strand="bottom",
        parameters=DesignParameters(max_tm=72),
    )
    product = results[0].valid_template_amplicons[0]
    assert product.start == 0
    assert product.end == len(TEMPLATE)


def test_evaluate_can_attach_bowtie_specificity(monkeypatch) -> None:
    def fake_assessment(outcomes, index, **options):
        assert index == "background"
        assert options["mismatches"] == 1
        outcome = outcomes[0]
        specificity = PrimerSpecificity("background", 1, 2, 3, 0, ())
        pairs = tuple(replace(pair, specificity=specificity) for pair in outcome.pairs)
        return [replace(outcome, pairs=pairs)]

    monkeypatch.setattr("smfprimer.evaluation.assess_specificity", fake_assessment)
    result = evaluate_primer_pairs(
        PrimerSet("known", FORWARD_ORDER, REVERSE_ORDER),
        TEMPLATE,
        workflow=Workflow.DEAMINASE,
        parameters=DesignParameters(max_tm=72),
        bowtie_index="background",
        mismatches=1,
    )[0]
    assert result.specificity is not None
    assert result.specificity.forward_hit_count == 2
    assert result.specificity.reverse_hit_count == 3
