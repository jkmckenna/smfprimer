"""Public data models used by the design API."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from typing import Any


class Workflow(StrEnum):
    DEAMINASE = "deaminase"
    CONVERSION = "conversion"


class ConvertedStrand(StrEnum):
    TOP = "top"
    BOTTOM = "bottom"


class TargetContext(StrEnum):
    CPG = "cpg"
    GPC = "gpc"
    BOTH = "both"


@dataclass(frozen=True)
class DesignParameters:
    """Primer3 candidate-generation and pair-ranking parameters."""

    min_length: int = 18
    optimum_length: int = 22
    max_length: int = 26
    min_tm: float = 50.0
    optimum_tm: float = 60.0
    max_tm: float = 68.0
    min_gc: float = 0.25
    max_gc: float = 0.75
    search_window: int = 120
    max_results: int = 5
    min_amplicon_size: int | None = None
    max_amplicon_size: int | None = None

    def __post_init__(self) -> None:
        if not 1 <= self.min_length <= self.optimum_length <= self.max_length:
            raise ValueError("primer lengths must satisfy 1 <= min <= optimum <= max")
        if not self.min_tm <= self.optimum_tm <= self.max_tm:
            raise ValueError("melting temperatures must satisfy min <= optimum <= max")
        if not 0 <= self.min_gc <= self.max_gc <= 1:
            raise ValueError("GC fractions must satisfy 0 <= min <= max <= 1")
        if self.search_window < self.max_length:
            raise ValueError("search_window must be at least max_length")
        if self.max_results < 1:
            raise ValueError("max_results must be positive")
        if self.min_amplicon_size is not None and self.min_amplicon_size < 1:
            raise ValueError("min_amplicon_size must be positive")
        if self.max_amplicon_size is not None and self.max_amplicon_size < 1:
            raise ValueError("max_amplicon_size must be positive")
        if (
            self.min_amplicon_size is not None
            and self.max_amplicon_size is not None
            and self.min_amplicon_size > self.max_amplicon_size
        ):
            raise ValueError("amplicon sizes must satisfy min <= max")


@dataclass(frozen=True)
class FeatureSegment:
    start: int
    end: int

    def __post_init__(self) -> None:
        if not 0 <= self.start < self.end:
            raise ValueError("feature segment must satisfy 0 <= start < end")


@dataclass(frozen=True)
class ReferenceFeature:
    """A coordinate-bearing annotation relative to a DesignTarget sequence."""

    type: str
    label: str
    segments: tuple[FeatureSegment, ...]
    strand: int = 0
    qualifiers: tuple[tuple[str, tuple[str, ...]], ...] = ()

    def __post_init__(self) -> None:
        if not self.type:
            raise ValueError("feature type must not be empty")
        if not self.segments:
            raise ValueError("feature must contain at least one segment")
        if self.strand not in {-1, 0, 1}:
            raise ValueError("feature strand must be -1, 0, or 1")

    @property
    def start(self) -> int:
        return min(segment.start for segment in self.segments)

    @property
    def end(self) -> int:
        return max(segment.end for segment in self.segments)

    @classmethod
    def create(
        cls,
        *,
        type: str,
        label: str = "",
        segments: list[tuple[int, int]] | tuple[tuple[int, int], ...],
        strand: int = 0,
        qualifiers: dict[str, Any] | None = None,
    ) -> ReferenceFeature:
        normalized = []
        for key, value in (qualifiers or {}).items():
            values = value if isinstance(value, (list, tuple)) else (value,)
            normalized.append((key, tuple(str(item) for item in values)))
        return cls(
            type=type,
            label=label,
            segments=tuple(FeatureSegment(start, end) for start, end in segments),
            strand=strand,
            qualifiers=tuple(sorted(normalized)),
        )

    def qualifier_dict(self) -> dict[str, tuple[str, ...]]:
        return dict(self.qualifiers)


@dataclass(frozen=True)
class DesignTarget:
    """A target interval within a supplied top-strand sequence window."""

    target_id: str
    reference_name: str
    sequence: str
    target_start: int
    target_end: int
    sequence_start: int = 0
    source: str = "sequence"
    feature_strand: str = "."
    metadata: tuple[tuple[str, str], ...] = ()
    reference_features: tuple[ReferenceFeature, ...] = ()
    topology: str = "linear"

    def __post_init__(self) -> None:
        if not self.target_id:
            raise ValueError("target_id must not be empty")
        if not 0 <= self.target_start < self.target_end <= len(self.sequence):
            raise ValueError("target must satisfy 0 <= start < end <= sequence length")
        if self.sequence_start < 0:
            raise ValueError("sequence_start must not be negative")
        if self.feature_strand not in {"+", "-", "."}:
            raise ValueError("feature_strand must be '+', '-', or '.'")
        if self.topology not in {"linear", "circular"}:
            raise ValueError("topology must be 'linear' or 'circular'")
        if any(feature.end > len(self.sequence) for feature in self.reference_features):
            raise ValueError("reference feature exceeds target sequence bounds")

    @property
    def reference_target_start(self) -> int:
        return self.sequence_start + self.target_start

    @property
    def reference_target_end(self) -> int:
        return self.sequence_start + self.target_end

    @classmethod
    def create(
        cls,
        *,
        target_id: str,
        reference_name: str,
        sequence: str,
        target_start: int,
        target_end: int,
        sequence_start: int = 0,
        source: str = "sequence",
        feature_strand: str = ".",
        metadata: dict[str, Any] | None = None,
        reference_features: tuple[ReferenceFeature, ...] = (),
        topology: str = "linear",
    ) -> DesignTarget:
        return cls(
            target_id=target_id,
            reference_name=reference_name,
            sequence=sequence,
            target_start=target_start,
            target_end=target_end,
            sequence_start=sequence_start,
            source=source,
            feature_strand=feature_strand,
            metadata=tuple(sorted((key, str(value)) for key, value in (metadata or {}).items())),
            reference_features=reference_features,
            topology=topology,
        )


@dataclass(frozen=True)
class CandidatePrimer:
    sequence: str
    unconverted_sequence: str
    side: str
    start: int
    end: int
    tm: float
    gc_fraction: float
    degeneracies: int
    score: float


@dataclass(frozen=True)
class OffTargetAmplicon:
    reference_name: str
    start: int
    end: int
    forward_strand: str
    forward_mismatches: int
    reverse_mismatches: int

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class PrimerSpecificity:
    index: str
    mismatches: int
    forward_hit_count: int
    reverse_hit_count: int
    intended_amplicon_count: int
    off_target_amplicons: tuple[OffTargetAmplicon, ...]

    @property
    def off_target_amplicon_count(self) -> int:
        return len(self.off_target_amplicons)

    @property
    def status(self) -> str:
        if self.off_target_amplicons:
            return "off_targets"
        if self.intended_amplicon_count:
            return "specific"
        return "no_off_targets"


@dataclass(frozen=True)
class PrimerPair:
    forward: CandidatePrimer
    reverse: CandidatePrimer
    amplicon_start: int
    amplicon_end: int
    converted_strand: ConvertedStrand
    score: float
    specificity: PrimerSpecificity | None = None

    @property
    def amplicon_length(self) -> int:
        return self.amplicon_end - self.amplicon_start


@dataclass(frozen=True)
class DesignOutcome:
    target: DesignTarget
    workflow: Workflow
    context: TargetContext
    converted_strand: ConvertedStrand
    parameters: DesignParameters
    pairs: tuple[PrimerPair, ...]
    message: str = ""

    @property
    def status(self) -> str:
        return "ok" if self.pairs else "no_candidates"
