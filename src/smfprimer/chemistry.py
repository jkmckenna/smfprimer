"""Sequence transformations for SMF primer chemistry."""

from __future__ import annotations

from .models import TargetContext, Workflow

_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def normalize_sequence(sequence: str) -> str:
    normalized = "".join(sequence.split()).upper()
    invalid = set(normalized) - set("ACGTRYSWKMBDHVN")
    if not normalized:
        raise ValueError("reference sequence is empty")
    if invalid:
        raise ValueError(f"reference contains unsupported bases: {''.join(sorted(invalid))}")
    return normalized


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def targeted_cytosines(sequence: str, context: TargetContext) -> frozenset[int]:
    """Return cytosines belonging to a requested context in a 5'-to-3' sequence."""
    cpg = {index for index in range(len(sequence) - 1) if sequence[index : index + 2] == "CG"}
    gpc = {index + 1 for index in range(len(sequence) - 1) if sequence[index : index + 2] == "GC"}
    if context is TargetContext.CPG:
        return frozenset(cpg)
    if context is TargetContext.GPC:
        return frozenset(gpc)
    return frozenset(cpg | gpc)


def encode_converted_template(
    sequence: str,
    workflow: Workflow,
    context: TargetContext = TargetContext.BOTH,
) -> str:
    """Encode one converted strand, using Y for an uncertain C/T position.

    In deaminase experiments every cytosine is treated as uncertain. In a
    conversion workflow, non-target cytosines have the expected converted T,
    while potentially methylated target-context cytosines use Y (C or T).
    """
    if workflow is Workflow.DEAMINASE:
        return sequence.replace("C", "Y")

    targets = targeted_cytosines(sequence, context)
    return "".join(
        "Y" if i in targets else "T" if base == "C" else base for i, base in enumerate(sequence)
    )
