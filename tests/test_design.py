from smfprimer import (
    ConvertedStrand,
    DesignParameters,
    Workflow,
    design_primers,
)

PARAMETERS = DesignParameters(
    min_length=4,
    optimum_length=4,
    max_length=4,
    min_tm=0,
    optimum_tm=12,
    max_tm=100,
    min_gc=0,
    max_gc=1,
    search_window=8,
    max_results=3,
)


def test_top_deaminase_forward_uses_y_and_reverse_uses_r() -> None:
    pairs = design_primers(
        "CCCCCCCCGGGGCCCC",
        8,
        10,
        workflow=Workflow.DEAMINASE,
        parameters=PARAMETERS,
    )
    assert pairs
    assert "Y" in pairs[0].forward.sequence
    assert "R" in pairs[0].reverse.sequence


def test_bottom_conversion_reports_coordinates_on_supplied_reference() -> None:
    pairs = design_primers(
        "AAAACCCCGGGGTTTTAAAA",
        8,
        12,
        converted_strand=ConvertedStrand.BOTTOM,
        parameters=PARAMETERS,
    )
    assert pairs
    assert pairs[0].amplicon_start <= 8
    assert pairs[0].amplicon_end >= 12


def test_bottom_deaminase_switches_reference_g_positions() -> None:
    pairs = design_primers(
        "GGGGGGGGAAAAGGGG",
        8,
        10,
        workflow=Workflow.DEAMINASE,
        converted_strand=ConvertedStrand.BOTTOM,
        parameters=PARAMETERS,
    )
    assert pairs
    assert "Y" in pairs[0].forward.sequence
    assert "R" in pairs[0].reverse.sequence


def test_degeneracy_is_preferred_over_other_scoring_terms() -> None:
    pairs = design_primers(
        "AAAACCCCAAAATTTTAAAAGGGG",
        12,
        16,
        workflow=Workflow.CONVERSION,
        parameters=PARAMETERS,
    )
    assert pairs[0].forward.degeneracies + pairs[0].reverse.degeneracies == 0


def test_amplicon_size_constraints_apply_to_complete_product() -> None:
    parameters = DesignParameters(
        min_length=4,
        optimum_length=4,
        max_length=4,
        min_tm=0,
        optimum_tm=12,
        max_tm=100,
        min_gc=0,
        max_gc=1,
        search_window=10,
        max_results=5,
        min_amplicon_size=14,
        max_amplicon_size=16,
    )
    pairs = design_primers("ACGT" * 10, 18, 22, parameters=parameters)
    assert pairs
    assert all(14 <= pair.amplicon_length <= 16 for pair in pairs)
