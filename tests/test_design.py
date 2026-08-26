from smfprimer import (
    ConvertedStrand,
    DesignParameters,
    Workflow,
    design_primers,
)
from smfprimer.chemistry import reverse_complement

FORWARD_SITE = "GAGAGATCTGGCAGCGGAGAG"
REVERSE_SITE = "CTTTTCTGTCACCAATCCTGTCCC"
TEMPLATE = FORWARD_SITE + "A" * 20 + reverse_complement(REVERSE_SITE)
PARAMETERS = DesignParameters(
    min_length=18,
    optimum_length=21,
    max_length=25,
    min_tm=35,
    optimum_tm=55,
    max_tm=75,
    min_gc=0.2,
    max_gc=0.8,
    search_window=30,
    max_results=3,
    min_amplicon_size=55,
    max_amplicon_size=70,
)


def test_top_deaminase_forward_uses_y_and_reverse_uses_r() -> None:
    pairs = design_primers(
        TEMPLATE,
        25,
        40,
        workflow=Workflow.DEAMINASE,
        parameters=PARAMETERS,
    )
    assert pairs
    assert "Y" in pairs[0].forward.sequence
    assert "R" in pairs[0].reverse.sequence


def test_bottom_conversion_reports_coordinates_on_supplied_reference() -> None:
    pairs = design_primers(
        reverse_complement(TEMPLATE),
        25,
        40,
        converted_strand=ConvertedStrand.BOTTOM,
        parameters=PARAMETERS,
    )
    assert pairs
    assert pairs[0].amplicon_start <= 25
    assert pairs[0].amplicon_end >= 40


def test_bottom_deaminase_switches_reference_g_positions() -> None:
    pairs = design_primers(
        reverse_complement(TEMPLATE),
        25,
        40,
        workflow=Workflow.DEAMINASE,
        converted_strand=ConvertedStrand.BOTTOM,
        parameters=PARAMETERS,
    )
    assert pairs
    assert "Y" in pairs[0].forward.sequence
    assert "R" in pairs[0].reverse.sequence


def test_primer3_penalty_replaces_custom_degeneracy_score() -> None:
    pairs = design_primers(
        TEMPLATE,
        25,
        40,
        workflow=Workflow.DEAMINASE,
        parameters=PARAMETERS,
    )
    assert pairs[0].forward.degeneracies + pairs[0].reverse.degeneracies > 0
    assert pairs[0].score < 100
    assert [pair.score for pair in pairs] == sorted(pair.score for pair in pairs)


def test_amplicon_size_constraints_apply_to_complete_product() -> None:
    parameters = DesignParameters(
        min_length=18,
        optimum_length=20,
        max_length=22,
        min_tm=30,
        optimum_tm=50,
        max_tm=80,
        min_gc=0,
        max_gc=1,
        search_window=30,
        max_results=5,
        min_amplicon_size=50,
        max_amplicon_size=60,
    )
    pairs = design_primers("ACGT" * 40, 60, 70, parameters=parameters)
    assert pairs
    assert all(50 <= pair.amplicon_length <= 60 for pair in pairs)
