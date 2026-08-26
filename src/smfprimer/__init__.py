"""Primer design for single-molecule footprinting experiments."""

from .design import design_primers, design_targets
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

__all__ = [
    "CandidatePrimer",
    "ConvertedStrand",
    "DesignOutcome",
    "DesignParameters",
    "DesignTarget",
    "PrimerPair",
    "TargetContext",
    "Workflow",
    "design_primers",
    "design_targets",
]

__version__ = "0.1.0"
