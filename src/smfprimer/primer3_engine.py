"""Small adapter around primer3-py's BoulderIO-style bindings."""

from __future__ import annotations

from typing import Any

from primer3 import bindings

from .chemistry import normalize_sequence
from .models import DesignParameters

_CONCRETE_PRIMER_BASE = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "A",
    "Y": "T",
    "S": "C",
    "W": "A",
    "K": "T",
    "M": "A",
    "B": "T",
    "D": "A",
    "H": "A",
    "V": "A",
    "N": "A",
}


def concrete_converted_template(sequence: str) -> str:
    """Return the expected converted allele accepted by Primer3.

    smfprimer-generated Y bases represent C/T conversion uncertainty and are
    represented by T for design. Other reference ambiguity remains unknown.
    """
    sequence = normalize_sequence(sequence)
    return "".join(base if base in "ACGT" else "T" if base == "Y" else "N" for base in sequence)


def concrete_primer(sequence: str) -> str:
    """Choose the lower-GC concrete member of an IUPAC oligo mixture."""
    return "".join(_CONCRETE_PRIMER_BASE[base] for base in normalize_sequence(sequence))


def global_args(parameters: DesignParameters, *, task: str = "generic") -> dict[str, Any]:
    """Translate smfprimer parameters into Primer3 global arguments."""
    return {
        "PRIMER_TASK": task,
        "PRIMER_PICK_LEFT_PRIMER": 1,
        "PRIMER_PICK_RIGHT_PRIMER": 1,
        "PRIMER_PICK_INTERNAL_OLIGO": 0,
        "PRIMER_NUM_RETURN": parameters.max_results,
        "PRIMER_MIN_SIZE": parameters.min_length,
        "PRIMER_OPT_SIZE": parameters.optimum_length,
        "PRIMER_MAX_SIZE": parameters.max_length,
        "PRIMER_MIN_TM": parameters.min_tm,
        "PRIMER_OPT_TM": parameters.optimum_tm,
        "PRIMER_MAX_TM": parameters.max_tm,
        "PRIMER_MIN_GC": parameters.min_gc * 100.0,
        "PRIMER_MAX_GC": parameters.max_gc * 100.0,
        "PRIMER_EXPLAIN_FLAG": 1,
    }


def design(
    template: str,
    target_start: int,
    target_end: int,
    parameters: DesignParameters,
) -> dict[str, Any]:
    """Run Primer3 pair design on a concrete converted template window."""
    minimum_product = max(
        parameters.min_amplicon_size or 1,
        (target_end - target_start) + 2 * parameters.min_length,
        parameters.max_length + 1,
    )
    maximum_product = min(parameters.max_amplicon_size or len(template), len(template))
    if minimum_product > maximum_product:
        return {"PRIMER_PAIR_NUM_RETURNED": 0}
    args = global_args(parameters)
    args["PRIMER_PRODUCT_SIZE_RANGE"] = [[minimum_product, maximum_product]]
    return bindings.design_primers(
        {
            "SEQUENCE_ID": "smfprimer_target",
            "SEQUENCE_TEMPLATE": concrete_converted_template(template),
            "SEQUENCE_TARGET": [target_start, target_end - target_start],
        },
        args,
    )


def check_pair(
    forward: str,
    reverse: str,
    parameters: DesignParameters,
) -> dict[str, Any]:
    """Have Primer3 calculate metrics and penalties for supplied primers."""
    args = global_args(parameters, task="check_primers")
    args["PRIMER_PICK_ANYWAY"] = 1
    args["PRIMER_NUM_RETURN"] = 1
    return bindings.design_primers(
        {
            "SEQUENCE_ID": "smfprimer_supplied_pair",
            "SEQUENCE_PRIMER": concrete_primer(forward),
            "SEQUENCE_PRIMER_REVCOMP": concrete_primer(reverse),
        },
        args,
    )


def problems(result: dict[str, Any], side: str) -> tuple[str, ...]:
    """Normalize Primer3's semicolon-delimited supplied-primer problems."""
    value = result.get(f"PRIMER_{side}_0_PROBLEMS", b"")
    if isinstance(value, bytes):
        value = value.decode()
    return tuple(part.strip() for part in str(value).split(";") if part.strip())
