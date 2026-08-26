from pathlib import Path

from smfprimer.cli import main

PERMISSIVE_OPTIONS = [
    "--product-size",
    "15:60",
    "--min-length",
    "4",
    "--optimum-length",
    "4",
    "--max-length",
    "4",
    "--min-tm",
    "0",
    "--max-tm",
    "100",
    "--min-gc",
    "0",
    "--max-gc",
    "1",
]


def test_fasta_cli_emits_one_result_for_each_record(tmp_path: Path, capsys) -> None:
    fasta = tmp_path / "targets.fa"
    fasta.write_text(">one\n" + "ACGT" * 20 + "\n>two\n" + "TGCA" * 20 + "\n")
    result = main(
        [
            "design",
            "fasta",
            str(fasta),
            *PERMISSIVE_OPTIONS,
        ]
    )
    lines = capsys.readouterr().out.splitlines()
    assert result == 0
    assert any(line.startswith("one\t") for line in lines[1:])
    assert any(line.startswith("two\t") for line in lines[1:])


def test_bed_cli_reports_genome_coordinates(tmp_path: Path, capsys) -> None:
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\n" + "ACGT" * 40 + "\n")
    bed = tmp_path / "targets.bed"
    bed.write_text("chr1\t60\t70\tbed-target\n")
    result = main(
        [
            "design",
            "bed",
            "--genome",
            str(genome),
            "--bed",
            str(bed),
            "--search-window",
            "20",
            *PERMISSIVE_OPTIONS,
        ]
    )
    output = capsys.readouterr().out
    assert result == 0
    assert "bed-target\tbed\tchr1" in output
    assert "\t60\t70\tok\t" in output


def test_tss_cli_accepts_gene_names(tmp_path: Path, capsys) -> None:
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\n" + "ACGT" * 60 + "\n")
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        "chr1\ttest\ttranscript\t101\t160\t.\t+\t.\t"
        'gene_id "gene1"; transcript_id "tx1"; gene_name "GENE";\n'
    )
    genes = tmp_path / "genes.txt"
    genes.write_text("GENE\n")
    result = main(
        [
            "design",
            "tss",
            "--genome",
            str(genome),
            "--gtf",
            str(gtf),
            "--genes",
            str(genes),
            "--tss-upstream",
            "2",
            "--tss-downstream",
            "2",
            "--search-window",
            "20",
            *PERMISSIVE_OPTIONS,
        ]
    )
    output = capsys.readouterr().out
    assert result == 0
    assert "GENE:tx1\ttss\tchr1" in output
