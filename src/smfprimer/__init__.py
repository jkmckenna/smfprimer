"""Primer design for single-molecule footprinting experiments."""

from .design import design_primers, design_targets
from .evaluation import (
    PrimerMetrics,
    PrimerPairEvaluation,
    PrimerSet,
    TemplateAmplicon,
    evaluate_primer_pairs,
)
from .genbank import format_genbank
from .models import (
    CandidatePrimer,
    ConvertedStrand,
    DesignOutcome,
    DesignParameters,
    DesignTarget,
    FeatureSegment,
    OffTargetAmplicon,
    PrimerPair,
    PrimerSpecificity,
    ReferenceFeature,
    TargetContext,
    Workflow,
)
from .specificity import assess_specificity

__all__ = [
    "CandidatePrimer",
    "ConvertedStrand",
    "DesignOutcome",
    "DesignParameters",
    "DesignTarget",
    "FeatureSegment",
    "OffTargetAmplicon",
    "PrimerPair",
    "PrimerMetrics",
    "PrimerPairEvaluation",
    "PrimerSet",
    "PrimerSpecificity",
    "ReferenceFeature",
    "TargetContext",
    "Workflow",
    "assess_specificity",
    "design_primers",
    "design_targets",
    "evaluate_primer_pairs",
    "format_genbank",
    "TemplateAmplicon",
]

__version__ = "0.2.0"
