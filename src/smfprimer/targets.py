"""Input adapters that normalize target sources for the design engine."""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

from .chemistry import normalize_sequence
from .models import DesignTarget


def read_fasta(path: str | Path) -> dict[str, str]:
    """Read a FASTA file, preserving record order and rejecting duplicate IDs."""
    records: dict[str, list[str]] = {}
    record_id: str | None = None
    with Path(path).open() as handle:
        for line_number, line in enumerate(handle, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                record_id = line[1:].split(maxsplit=1)[0]
                if not record_id:
                    raise ValueError(f"empty FASTA identifier at line {line_number}")
                if record_id in records:
                    raise ValueError(f"duplicate FASTA identifier: {record_id}")
                records[record_id] = []
            elif record_id is None:
                raise ValueError(f"FASTA sequence before first header at line {line_number}")
            else:
                records[record_id].append(line)
    if not records:
        raise ValueError("FASTA file contains no records")
    return {record_id: normalize_sequence("".join(lines)) for record_id, lines in records.items()}


@dataclass(frozen=True)
class _FastaIndexEntry:
    length: int
    offset: int
    line_bases: int
    line_width: int


class IndexedFasta:
    """Small read-only FASTA index supporting genome-scale random access.

    The index is built in memory without loading chromosome sequences or
    writing files beside the user's reference.
    """

    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)
        self._index = self._build_index()

    @property
    def references(self) -> tuple[str, ...]:
        return tuple(self._index)

    def length(self, reference: str) -> int:
        try:
            return self._index[reference].length
        except KeyError as error:
            raise ValueError(f"reference {reference!r} is absent from {self.path}") from error

    def fetch(self, reference: str, start: int, end: int) -> str:
        try:
            entry = self._index[reference]
        except KeyError as error:
            raise ValueError(f"reference {reference!r} is absent from {self.path}") from error
        if not 0 <= start <= end <= entry.length:
            raise ValueError(
                f"invalid interval {reference}:{start}-{end}; reference length is {entry.length}"
            )
        chunks: list[bytes] = []
        position = start
        with self.path.open("rb") as handle:
            while position < end:
                line_index, within_line = divmod(position, entry.line_bases)
                count = min(end - position, entry.line_bases - within_line)
                handle.seek(entry.offset + line_index * entry.line_width + within_line)
                chunk = handle.read(count)
                if len(chunk) != count:
                    raise ValueError(f"unexpected end of FASTA while reading {reference}")
                chunks.append(chunk)
                position += count
        return normalize_sequence(b"".join(chunks).decode("ascii"))

    def _build_index(self) -> dict[str, _FastaIndexEntry]:
        index: dict[str, _FastaIndexEntry] = {}
        current_id: str | None = None
        offset = length = line_bases = line_width = previous_bases = 0

        def finish() -> None:
            if current_id is None:
                return
            if length == 0:
                raise ValueError(f"FASTA record {current_id!r} is empty")
            index[current_id] = _FastaIndexEntry(length, offset, line_bases, line_width)

        with self.path.open("rb") as handle:
            while line := handle.readline():
                if line.startswith(b">"):
                    finish()
                    current_id = line[1:].decode("utf-8").strip().split(maxsplit=1)[0]
                    if not current_id:
                        raise ValueError("empty FASTA identifier")
                    if current_id in index:
                        raise ValueError(f"duplicate FASTA identifier: {current_id}")
                    offset = handle.tell()
                    length = line_bases = line_width = previous_bases = 0
                    continue
                if current_id is None:
                    if line.strip():
                        raise ValueError("FASTA sequence appears before its first header")
                    continue
                stripped = line.rstrip(b"\r\n")
                if not stripped:
                    raise ValueError("blank lines within FASTA records are not supported")
                if line_bases == 0:
                    line_bases, line_width = len(stripped), len(line)
                elif previous_bases != line_bases:
                    raise ValueError("only the final sequence line of a FASTA record may be short")
                elif len(stripped) > line_bases or (
                    len(stripped) == line_bases and len(line) not in {line_bases, line_width}
                ):
                    raise ValueError("FASTA sequence lines must use a consistent width")
                previous_bases = len(stripped)
                length += previous_bases
        finish()
        if not index:
            raise ValueError("FASTA file contains no records")
        return index


def sequence_target(
    sequence: str,
    target_start: int,
    target_end: int,
    *,
    target_id: str = "target",
) -> DesignTarget:
    return DesignTarget.create(
        target_id=target_id,
        reference_name=target_id,
        sequence=normalize_sequence(sequence),
        target_start=target_start,
        target_end=target_end,
    )


def fasta_targets(path: str | Path, *, target_width: int = 1) -> list[DesignTarget]:
    if target_width < 1:
        raise ValueError("target_width must be positive")
    targets = []
    for record_id, sequence in read_fasta(path).items():
        if target_width > len(sequence):
            raise ValueError(f"target_width exceeds FASTA record length for {record_id!r}")
        start = (len(sequence) - target_width) // 2
        targets.append(
            DesignTarget.create(
                target_id=record_id,
                reference_name=record_id,
                sequence=sequence,
                target_start=start,
                target_end=start + target_width,
                source="fasta",
            )
        )
    return targets


def bed_targets(
    genome_path: str | Path,
    bed_path: str | Path,
    *,
    flank: int,
) -> list[DesignTarget]:
    if flank < 1:
        raise ValueError("flank must be positive")
    genome = IndexedFasta(genome_path)
    targets = []
    with Path(bed_path).open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"BED line {line_number} has fewer than three columns")
            reference, start_text, end_text = fields[:3]
            try:
                start, end = int(start_text), int(end_text)
            except ValueError as error:
                raise ValueError(f"BED line {line_number} has non-integer coordinates") from error
            reference_length = genome.length(reference)
            if not 0 <= start < end <= reference_length:
                raise ValueError(f"BED line {line_number} is outside {reference}")
            window_start = max(0, start - flank)
            window_end = min(reference_length, end + flank)
            target_id = (
                fields[3] if len(fields) >= 4 and fields[3] else f"{reference}:{start}-{end}"
            )
            strand = fields[5] if len(fields) >= 6 and fields[5] in {"+", "-"} else "."
            targets.append(
                DesignTarget.create(
                    target_id=target_id,
                    reference_name=reference,
                    sequence=genome.fetch(reference, window_start, window_end),
                    target_start=start - window_start,
                    target_end=end - window_start,
                    sequence_start=window_start,
                    source="bed",
                    feature_strand=strand,
                )
            )
    if not targets:
        raise ValueError("BED file contains no target intervals")
    return targets


@dataclass(frozen=True)
class _Transcript:
    gene_id: str
    gene_name: str
    transcript_id: str
    reference: str
    start: int
    end: int
    strand: str

    @property
    def tss(self) -> int:
        return self.start if self.strand == "+" else self.end - 1


def _gtf_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in text.split(";"):
        item = item.strip()
        if not item:
            continue
        key, separator, value = item.partition(" ")
        if not separator:
            key, separator, value = item.partition("=")
        if separator:
            attributes[key] = value.strip().strip('"')
    return attributes


def _transcripts(gtf_path: str | Path) -> Iterator[_Transcript]:
    with Path(gtf_path).open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"GTF line {line_number} does not have nine columns")
            reference, _, feature, start_text, end_text, _, strand, _, attribute_text = fields
            if feature != "transcript":
                continue
            if strand not in {"+", "-"}:
                continue
            attributes = _gtf_attributes(attribute_text)
            gene_id = attributes.get("gene_id", "")
            transcript_id = attributes.get("transcript_id", "")
            if not gene_id or not transcript_id:
                continue
            yield _Transcript(
                gene_id=gene_id,
                gene_name=attributes.get("gene_name", gene_id),
                transcript_id=transcript_id,
                reference=reference,
                start=int(start_text) - 1,
                end=int(end_text),
                strand=strand,
            )


def _gene_names(path: str | Path) -> list[str]:
    names = [
        line.strip()
        for line in Path(path).read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]
    if not names:
        raise ValueError("gene list contains no identifiers")
    return names


def tss_targets(
    genome_path: str | Path,
    gtf_path: str | Path,
    genes_path: str | Path,
    *,
    upstream: int,
    downstream: int,
    flank: int,
    policy: str = "all",
) -> list[DesignTarget]:
    if min(upstream, downstream) < 0:
        raise ValueError("TSS upstream/downstream distances must not be negative")
    if policy not in {"all", "longest"}:
        raise ValueError("TSS policy must be 'all' or 'longest'")
    requested = _gene_names(genes_path)
    requested_set = set(requested)
    matches: dict[str, list[_Transcript]] = {name: [] for name in requested}
    for transcript in _transcripts(gtf_path):
        for identifier in {transcript.gene_id, transcript.gene_name} & requested_set:
            matches[identifier].append(transcript)
    missing = [name for name in requested if not matches[name]]
    if missing:
        raise ValueError("genes absent from transcript features in GTF: " + ", ".join(missing))

    genome = IndexedFasta(genome_path)
    targets: list[DesignTarget] = []
    for requested_name in requested:
        transcripts = matches[requested_name]
        if policy == "longest":
            transcripts = [max(transcripts, key=lambda item: item.end - item.start)]
        else:
            unique: dict[tuple[str, int, str], _Transcript] = {}
            for transcript in transcripts:
                unique.setdefault(
                    (transcript.reference, transcript.tss, transcript.strand), transcript
                )
            transcripts = list(unique.values())

        for transcript in transcripts:
            if transcript.strand == "+":
                target_start = transcript.tss - upstream
                target_end = transcript.tss + downstream + 1
            else:
                target_start = transcript.tss - downstream
                target_end = transcript.tss + upstream + 1
            reference_length = genome.length(transcript.reference)
            if target_start < 0 or target_end > reference_length:
                raise ValueError(f"TSS window for {requested_name!r} exceeds chromosome bounds")
            window_start = max(0, target_start - flank)
            window_end = min(reference_length, target_end + flank)
            targets.append(
                DesignTarget.create(
                    target_id=f"{requested_name}:{transcript.transcript_id}",
                    reference_name=transcript.reference,
                    sequence=genome.fetch(transcript.reference, window_start, window_end),
                    target_start=target_start - window_start,
                    target_end=target_end - window_start,
                    sequence_start=window_start,
                    source="tss",
                    feature_strand=transcript.strand,
                    metadata={
                        "gene_id": transcript.gene_id,
                        "gene_name": transcript.gene_name,
                        "transcript_id": transcript.transcript_id,
                        "tss": transcript.tss,
                    },
                )
            )
    return targets
