import pytest

from smfprimer.chemistry import encode_converted_template, reverse_complement, targeted_cytosines
from smfprimer.models import TargetContext, Workflow


def test_deaminase_marks_every_cytosine_as_c_or_t() -> None:
    assert encode_converted_template("ACGCC", Workflow.DEAMINASE) == "AYGYY"


@pytest.mark.parametrize(
    ("context", "expected"),
    [
        (TargetContext.CPG, frozenset({1})),
        (TargetContext.GPC, frozenset({4})),
        (TargetContext.BOTH, frozenset({1, 4})),
    ],
)
def test_target_contexts(context: TargetContext, expected: frozenset[int]) -> None:
    assert targeted_cytosines("ACGGCAT", context) == expected


def test_conversion_changes_non_target_cytosines_and_marks_targets() -> None:
    assert encode_converted_template("ACGGCAC", Workflow.CONVERSION, TargetContext.CPG) == "AYGGTAT"


def test_reverse_complement_preserves_iupac_ambiguity() -> None:
    assert reverse_complement("ATYGR") == "YCRAT"
