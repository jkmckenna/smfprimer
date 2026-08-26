from dataclasses import replace
from io import StringIO

from Bio import SeqIO

from smfprimer import DesignParameters, PrimerSpecificity, design_targets, format_genbank
from smfprimer.models import ReferenceFeature
from smfprimer.targets import sequence_target

PARAMETERS = DesignParameters(
    min_length=18,
    optimum_length=20,
    max_length=22,
    min_tm=30,
    optimum_tm=50,
    max_tm=80,
    min_gc=0,
    max_gc=1,
    search_window=30,
    max_results=1,
    min_amplicon_size=50,
    max_amplicon_size=70,
)


def test_genbank_contains_ranked_amplicon_and_primer_features() -> None:
    target = sequence_target("ACGT" * 40, 60, 70, target_id="locus one")
    outcomes = design_targets([target], parameters=PARAMETERS)
    pair = outcomes[0].pairs[0]
    content = format_genbank(outcomes)

    assert content.startswith("LOCUS       locus_one_1")
    assert '                     /label="locus one_required_target"' in content
    assert '                     /label="rank_1_amplicon"' in content
    assert '                     /label="rank_1_forward_primer"' in content
    assert '                     /label="rank_1_reverse_primer"' in content
    assert f"     misc_feature    {pair.amplicon_start + 1}..{pair.amplicon_end}" in content
    assert f"     primer_bind     {pair.forward.start + 1}..{pair.forward.end}" in content
    assert (
        f"     primer_bind     complement({pair.reverse.start + 1}..{pair.reverse.end})" in content
    )
    assert "ORIGIN\n        1 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt" in content


def test_bottom_strand_reverses_primer_feature_arrows() -> None:
    target = sequence_target("ACGT" * 40, 60, 70, target_id="bottom")
    outcomes = design_targets(
        [target],
        converted_strand="bottom",
        parameters=PARAMETERS,
    )
    pair = outcomes[0].pairs[0]
    content = format_genbank(outcomes)
    assert (
        f"     primer_bind     complement({pair.forward.start + 1}..{pair.forward.end})" in content
    )
    assert f"     primer_bind     {pair.reverse.start + 1}..{pair.reverse.end}" in content


def test_genbank_writes_one_record_per_outcome() -> None:
    targets = [
        sequence_target("ACGT" * 40, 60, 70, target_id="one"),
        sequence_target("TGCA" * 40, 60, 70, target_id="two"),
    ]
    content = format_genbank(design_targets(targets, parameters=PARAMETERS))
    assert content.count("LOCUS       ") == 2
    assert content.count("\n//\n") == 2


def test_genbank_does_not_split_path_qualifiers() -> None:
    target = sequence_target("ACGT" * 40, 60, 70, target_id="locus")
    outcome = design_targets([target], parameters=PARAMETERS)[0]
    index = "/long/path/without/spaces/to/a/GRCh38_bowtie1_index"
    specificity = PrimerSpecificity(index, 2, 1, 1, 1, ())
    pair = replace(outcome.pairs[0], specificity=specificity)
    content = format_genbank([replace(outcome, pairs=(pair,))])
    assert f'/specificity_index="{index}"' in content


def test_genbank_preserves_long_source_feature_labels_exactly() -> None:
    label = "feature prefix " + "long_identifier_" * 8
    feature = ReferenceFeature.create(
        type="promoter",
        label=label,
        segments=[(2, 8)],
        qualifiers={"label": label},
    )
    target = sequence_target("ACGT" * 40, 60, 70, target_id="locus")
    target = replace(target, reference_features=(feature,))
    outcome = design_targets([target], parameters=PARAMETERS)[0]
    parsed = SeqIO.read(StringIO(format_genbank([outcome])), "genbank")
    promoter = next(item for item in parsed.features if item.type == "promoter")
    assert promoter.qualifiers["label"] == [label]
