from pathlib import Path

import pytest

from smfprimer.targets import IndexedFasta, bed_targets, fasta_targets, tss_targets


def _write(path: Path, content: str) -> Path:
    path.write_text(content)
    return path


def test_indexed_fasta_fetches_wrapped_sequence_without_loading_it(tmp_path: Path) -> None:
    fasta = _write(tmp_path / "genome.fa", ">chr1 description\nACGTN\nTGCAA\nCC\n")
    genome = IndexedFasta(fasta)
    assert genome.references == ("chr1",)
    assert genome.length("chr1") == 12
    assert genome.fetch("chr1", 3, 10) == "TNTGCAA"


def test_multi_fasta_creates_one_centered_target_per_record(tmp_path: Path) -> None:
    fasta = _write(tmp_path / "targets.fa", ">one\nAACCGGTT\n>two\nTTTTCCCCAA\n")
    targets = fasta_targets(fasta, target_width=2)
    assert [target.target_id for target in targets] == ["one", "two"]
    assert [(target.target_start, target.target_end) for target in targets] == [(3, 5), (4, 6)]


def test_bed_target_preserves_required_span_and_reference_coordinates(tmp_path: Path) -> None:
    genome = _write(tmp_path / "genome.fa", ">chr1\n" + "ACGT" * 20 + "\n")
    bed = _write(tmp_path / "targets.bed", "chr1\t20\t30\tlocus-a\t0\t-\n")
    target = bed_targets(genome, bed, flank=8)[0]
    assert target.target_id == "locus-a"
    assert target.feature_strand == "-"
    assert target.sequence_start == 12
    assert (target.target_start, target.target_end) == (8, 18)
    assert (target.reference_target_start, target.reference_target_end) == (20, 30)


def test_tss_targets_are_strand_aware_and_match_gene_id_or_name(tmp_path: Path) -> None:
    genome = _write(tmp_path / "genome.fa", ">chr1\n" + "ACGT" * 60 + "\n")
    gtf = _write(
        tmp_path / "genes.gtf",
        "chr1\ttest\ttranscript\t51\t100\t.\t+\t.\t"
        'gene_id "gene1"; transcript_id "tx1"; gene_name "PLUS";\n'
        "chr1\ttest\ttranscript\t121\t160\t.\t-\t.\t"
        'gene_id "gene2"; transcript_id "tx2"; gene_name "MINUS";\n',
    )
    genes = _write(tmp_path / "genes.txt", "gene1\nMINUS\n")
    targets = tss_targets(
        genome,
        gtf,
        genes,
        upstream=2,
        downstream=3,
        flank=10,
    )
    assert [(target.reference_target_start, target.reference_target_end) for target in targets] == [
        (48, 54),
        (156, 162),
    ]
    assert [target.feature_strand for target in targets] == ["+", "-"]


def test_tss_loader_reports_missing_genes(tmp_path: Path) -> None:
    genome = _write(tmp_path / "genome.fa", ">chr1\n" + "A" * 100 + "\n")
    gtf = _write(tmp_path / "genes.gtf", "")
    genes = _write(tmp_path / "genes.txt", "missing\n")
    with pytest.raises(ValueError, match="missing"):
        tss_targets(genome, gtf, genes, upstream=1, downstream=1, flank=5)
