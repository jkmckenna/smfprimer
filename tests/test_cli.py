from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from smfprimer.cli import main

PERMISSIVE_OPTIONS = [
    "--product-size",
    "45:80",
    "--min-length",
    "18",
    "--optimum-length",
    "20",
    "--max-length",
    "22",
    "--min-tm",
    "30",
    "--optimum-tm",
    "50",
    "--max-tm",
    "80",
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
            "25",
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
            "25",
            *PERMISSIVE_OPTIONS,
        ]
    )
    output = capsys.readouterr().out
    assert result == 0
    assert "GENE:tx1\ttss\tchr1" in output


def test_cli_forwards_bowtie_specificity_options(tmp_path: Path, capsys, monkeypatch) -> None:
    fasta = tmp_path / "target.fa"
    fasta.write_text(">one\n" + "ACGT" * 20 + "\n")
    received = {}

    def fake_assessment(outcomes, index, **options):
        received["index"] = index
        received.update(options)
        return outcomes

    monkeypatch.setattr("smfprimer.cli.assess_specificity", fake_assessment)
    result = main(
        [
            "design",
            "fasta",
            str(fasta),
            *PERMISSIVE_OPTIONS,
            "--bowtie-index",
            "genome-index",
            "--bowtie-executable",
            "/bin/bowtie",
            "--specificity-mismatches",
            "1",
            "--specificity-max-expansions",
            "64",
        ]
    )
    capsys.readouterr()
    assert result == 0
    assert received == {
        "index": "genome-index",
        "mismatches": 1,
        "bowtie": "/bin/bowtie",
        "maximum_expansions": 64,
    }


def test_file_output_gets_automatic_genbank_companion(tmp_path: Path) -> None:
    fasta = tmp_path / "target.fa"
    fasta.write_text(">one\n" + "ACGT" * 20 + "\n")
    output = tmp_path / "primers.tsv"
    result = main(
        [
            "design",
            "fasta",
            str(fasta),
            *PERMISSIVE_OPTIONS,
            "--output",
            str(output),
        ]
    )
    assert result == 0
    assert output.is_file()
    genbank = tmp_path / "primers.gb"
    assert genbank.is_file()
    assert '                     /label="rank_1_amplicon"' in genbank.read_text()


def test_no_genbank_disables_automatic_companion(tmp_path: Path) -> None:
    fasta = tmp_path / "target.fa"
    fasta.write_text(">one\n" + "ACGT" * 20 + "\n")
    output = tmp_path / "primers.tsv"
    result = main(
        [
            "design",
            "fasta",
            str(fasta),
            *PERMISSIVE_OPTIONS,
            "--output",
            str(output),
            "--no-genbank",
        ]
    )
    assert result == 0
    assert output.is_file()
    assert not (tmp_path / "primers.gb").exists()


def test_annotated_cli_uses_configured_feature_and_preserves_labels(tmp_path: Path) -> None:
    record = SeqRecord(Seq("ACGT" * 30), id="plasmid", name="plasmid")
    record.annotations = {"molecule_type": "DNA", "topology": "circular"}
    record.features = [
        SeqFeature(
            SimpleLocation(30, 40),
            type="misc_feature",
            qualifiers={"label": ["capture_here"]},
        ),
        SeqFeature(
            SimpleLocation(10, 20),
            type="promoter",
            qualifiers={"label": ["original_promoter"]},
        ),
    ]
    annotated = tmp_path / "input.gb"
    SeqIO.write(record, annotated, "genbank")
    names = tmp_path / "required_features.txt"
    names.write_text("capture_here\n")
    output = tmp_path / "primers.tsv"
    result = main(
        [
            "design",
            "annotated",
            str(annotated),
            "--required-features-file",
            str(names),
            *PERMISSIVE_OPTIONS,
            "--output",
            str(output),
        ]
    )
    assert result == 0
    records = list(SeqIO.parse(tmp_path / "primers.gb", "genbank"))
    assert len(records) == 1
    labels = {feature.qualifiers.get("label", [""])[0] for feature in records[0].features}
    assert "capture_here" in labels
    assert "original_promoter" in labels
    assert any(label.endswith("_required_target") for label in labels)
    assert "rank_1_forward_primer" in labels
