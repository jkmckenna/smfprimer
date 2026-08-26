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
    """Candidate-generation and ranking parameters.

    Melting temperatures use a dependency-free Wallace-rule estimate. They
    are intended for candidate screening rather than final assay validation.
    """

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

    def __post_init__(self) -> None:
        if not self.target_id:
            raise ValueError("target_id must not be empty")
        if not 0 <= self.target_start < self.target_end <= len(self.sequence):
            raise ValueError("target must satisfy 0 <= start < end <= sequence length")
        if self.sequence_start < 0:
            raise ValueError("sequence_start must not be negative")
        if self.feature_strand not in {"+", "-", "."}:
            raise ValueError("feature_strand must be '+', '-', or '.'")

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
class PrimerPair:
    forward: CandidatePrimer
    reverse: CandidatePrimer
    amplicon_start: int
    amplicon_end: int
    converted_strand: ConvertedStrand
    score: float

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
