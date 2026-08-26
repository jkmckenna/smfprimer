"""Stable serialization of primer design outcomes."""

from __future__ import annotations

import csv
import json
from dataclasses import asdict
from io import StringIO

from .chemistry import encode_converted_template, reverse_complement
from .models import ConvertedStrand, DesignOutcome, PrimerPair

FIELDS = [
    "target_id",
    "source",
    "reference_name",
    "feature_strand",
    "target_start",
    "target_end",
    "status",
    "message",
    "rank",
    "workflow",
    "context",
    "converted_strand",
    "amplicon_start",
    "amplicon_end",
    "amplicon_length",
    "amplicon_reference_sequence",
    "amplicon_unconverted_sequence",
    "amplicon_converted_sequence",
    "forward_start",
    "forward_end",
    "forward_unconverted_sequence",
    "forward_order_sequence",
    "forward_tm",
    "forward_gc_fraction",
    "forward_degeneracies",
    "reverse_start",
    "reverse_end",
    "reverse_unconverted_sequence",
    "reverse_order_sequence",
    "reverse_tm",
    "reverse_gc_fraction",
    "reverse_degeneracies",
    "pair_score",
    "metadata",
    "parameters",
]


def outcome_rows(outcomes: list[DesignOutcome]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for outcome in outcomes:
        if not outcome.pairs:
            rows.append(_row(outcome, None, None))
            continue
        rows.extend(_row(outcome, pair, rank) for rank, pair in enumerate(outcome.pairs, 1))
    return rows


def _row(
    outcome: DesignOutcome,
    pair: PrimerPair | None,
    rank: int | None,
) -> dict[str, object]:
    target = outcome.target
    row: dict[str, object] = {field: None for field in FIELDS}
    row.update(
        {
            "target_id": target.target_id,
            "source": target.source,
            "reference_name": target.reference_name,
            "feature_strand": target.feature_strand,
            "target_start": target.reference_target_start,
            "target_end": target.reference_target_end,
            "status": outcome.status,
            "message": outcome.message,
            "rank": rank,
            "workflow": outcome.workflow.value,
            "context": outcome.context.value,
            "converted_strand": outcome.converted_strand.value,
            "metadata": json.dumps(dict(target.metadata), sort_keys=True),
            "parameters": json.dumps(asdict(outcome.parameters), sort_keys=True),
        }
    )
    if pair is None:
        return row

    local_start, local_end = pair.amplicon_start, pair.amplicon_end
    reference_amplicon = target.sequence[local_start:local_end]
    if outcome.converted_strand is ConvertedStrand.TOP:
        unconverted_amplicon = reference_amplicon
    else:
        unconverted_amplicon = reverse_complement(reference_amplicon)
    converted_amplicon = encode_converted_template(
        unconverted_amplicon, outcome.workflow, outcome.context
    )
    row.update(
        {
            "amplicon_start": target.sequence_start + local_start,
            "amplicon_end": target.sequence_start + local_end,
            "amplicon_length": pair.amplicon_length,
            "amplicon_reference_sequence": reference_amplicon,
            "amplicon_unconverted_sequence": unconverted_amplicon,
            "amplicon_converted_sequence": converted_amplicon,
            "forward_start": target.sequence_start + pair.forward.start,
            "forward_end": target.sequence_start + pair.forward.end,
            "forward_unconverted_sequence": pair.forward.unconverted_sequence,
            "forward_order_sequence": pair.forward.sequence,
            "forward_tm": round(pair.forward.tm, 3),
            "forward_gc_fraction": round(pair.forward.gc_fraction, 5),
            "forward_degeneracies": pair.forward.degeneracies,
            "reverse_start": target.sequence_start + pair.reverse.start,
            "reverse_end": target.sequence_start + pair.reverse.end,
            "reverse_unconverted_sequence": pair.reverse.unconverted_sequence,
            "reverse_order_sequence": pair.reverse.sequence,
            "reverse_tm": round(pair.reverse.tm, 3),
            "reverse_gc_fraction": round(pair.reverse.gc_fraction, 5),
            "reverse_degeneracies": pair.reverse.degeneracies,
            "pair_score": round(pair.score, 3),
        }
    )
    return row


def format_tsv(outcomes: list[DesignOutcome]) -> str:
    stream = StringIO()
    writer = csv.DictWriter(stream, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(outcome_rows(outcomes))
    return stream.getvalue()


def format_json(outcomes: list[DesignOutcome]) -> str:
    rows = outcome_rows(outcomes)
    for row in rows:
        row["metadata"] = json.loads(str(row["metadata"]))
        row["parameters"] = json.loads(str(row["parameters"]))
    return json.dumps(rows, indent=2) + "\n"
