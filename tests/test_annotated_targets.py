import gzip
import struct
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from smfprimer.targets import (
    annotate_targets_from_gtf,
    annotated_targets,
    bed_targets,
)


def _genbank(path: Path) -> Path:
    record = SeqRecord(Seq("ACGT" * 30), id="plasmid", name="plasmid")
    record.annotations = {"molecule_type": "DNA", "topology": "circular"}
    record.features = [
        SeqFeature(
            SimpleLocation(20, 30, strand=1),
            type="misc_feature",
            qualifiers={"label": ["required_interval"]},
        ),
        SeqFeature(
            SimpleLocation(60, 72, strand=-1),
            type="misc_feature",
            qualifiers={"label": ["required_interval"]},
        ),
        SeqFeature(
            SimpleLocation(35, 45, strand=1),
            type="promoter",
            qualifiers={"label": ["retained_promoter"], "note": ["keep this label"]},
        ),
    ]
    SeqIO.write(record, path, "genbank")
    return path


def _packet(packet_type: int, data: bytes) -> bytes:
    return bytes([packet_type]) + struct.pack(">I", len(data)) + data


def _snapgene(path: Path) -> Path:
    sequence = ("ACGT" * 30).encode()
    cookie = struct.pack(">8sHHH", b"SnapGene", 1, 15, 17)
    features = (
        b'<?xml version="1.0"?><Features>'
        b'<Feature name="required_interval" type="misc_feature" directionality="1">'
        b'<Segment range="21-30" type="standard"/>'
        b"</Feature>"
        b'<Feature name="retained_promoter" type="promoter" directionality="1">'
        b'<Segment range="36-45" type="standard"/>'
        b"</Feature>"
        b"</Features>"
    )
    path.write_bytes(
        _packet(0x09, cookie) + _packet(0x00, bytes([1]) + sequence) + _packet(0x0A, features)
    )
    return path


def test_genbank_input_creates_one_target_per_required_feature(tmp_path: Path) -> None:
    targets = annotated_targets(_genbank(tmp_path / "input.gb"))
    assert [(target.target_start, target.target_end) for target in targets] == [
        (20, 30),
        (60, 72),
    ]
    assert [target.feature_strand for target in targets] == ["+", "-"]
    assert all(target.topology == "circular" for target in targets)
    assert all(
        "retained_promoter" in {feature.label for feature in target.reference_features}
        for target in targets
    )


def test_snapgene_input_preserves_features(tmp_path: Path) -> None:
    targets = annotated_targets(_snapgene(tmp_path / "input.dna"))
    assert len(targets) == 1
    assert (targets[0].target_start, targets[0].target_end) == (20, 30)
    assert targets[0].topology == "circular"
    assert {feature.label for feature in targets[0].reference_features} == {
        "required_interval",
        "retained_promoter",
    }


def test_configured_feature_names_replace_the_default(tmp_path: Path) -> None:
    path = _genbank(tmp_path / "input.gb")
    with pytest.raises(ValueError, match="no required features"):
        annotated_targets(path, required_feature_names=["capture_A"])
    targets = annotated_targets(
        path,
        required_feature_names=["retained_promoter"],
    )
    assert len(targets) == 1
    assert (targets[0].target_start, targets[0].target_end) == (35, 45)


def test_gtf_gz_features_are_clipped_and_remapped_to_bed_window(tmp_path: Path) -> None:
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\n" + "ACGT" * 50 + "\n")
    bed = tmp_path / "targets.bed"
    bed.write_text("chr1\t50\t60\ttarget\n")
    targets = bed_targets(genome, bed, flank=10)
    gtf = tmp_path / "features.gtf.gz"
    with gzip.open(gtf, "wt") as handle:
        handle.write('chr1\ttest\tgene\t36\t55\t.\t+\t.\tgene_id "gene1"; gene_name "CLIPPED";\n')
        handle.write('chr1\ttest\texon\t56\t65\t.\t-\t.\tgene_id "gene1"; exon_id "exon1";\n')
        handle.write(
            'chr1\ttest\tgene\t101\t120\t.\t+\t.\tgene_id "outside"; gene_name "OUTSIDE";\n'
        )
    annotated = annotate_targets_from_gtf(targets, gtf)[0]
    assert [
        (feature.type, feature.label, feature.start, feature.end, feature.strand)
        for feature in annotated.reference_features
    ] == [
        ("gene", "CLIPPED", 0, 15, 1),
        ("exon", "exon1", 15, 25, -1),
    ]
