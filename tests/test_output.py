import json

from smfprimer import DesignParameters, Workflow, design_targets
from smfprimer.chemistry import reverse_complement
from smfprimer.output import format_json, format_tsv, outcome_rows
from smfprimer.targets import sequence_target

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
    max_results=1,
    min_amplicon_size=55,
    max_amplicon_size=70,
)

FORWARD_SITE = "GAGAGATCTGGCAGCGGAGAG"
REVERSE_SITE = "CTTTTCTGTCACCAATCCTGTCCC"
TEMPLATE = FORWARD_SITE + "A" * 20 + reverse_complement(REVERSE_SITE)


def test_output_distinguishes_reference_unconverted_and_order_sequences() -> None:
    target = sequence_target(TEMPLATE, 25, 40, target_id="locus")
    outcomes = design_targets([target], workflow=Workflow.DEAMINASE, parameters=PARAMETERS)
    row = outcome_rows(outcomes)[0]
    assert row["status"] == "ok"
    assert row["amplicon_reference_sequence"]
    assert row["amplicon_unconverted_sequence"]
    assert "Y" in row["forward_order_sequence"]
    assert "Y" not in row["forward_unconverted_sequence"]
    assert "R" in row["reverse_order_sequence"]
    assert "R" not in row["reverse_unconverted_sequence"]
    assert format_tsv(outcomes).startswith("target_id\tsource\t")
    assert json.loads(format_json(outcomes))[0]["target_id"] == "locus"


def test_failed_target_is_retained_in_batch_output() -> None:
    target = sequence_target("A" * 30, 14, 16)
    parameters = DesignParameters(min_gc=0.9, max_gc=1.0)
    row = outcome_rows(design_targets([target], parameters=parameters))[0]
    assert row["status"] == "no_candidates"
    assert "no primer pairs" in row["message"]
