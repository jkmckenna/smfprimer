import json
from pathlib import Path

from smfprimer import (
    CandidatePrimer,
    ConvertedStrand,
    DesignOutcome,
    DesignParameters,
    DesignTarget,
    PrimerPair,
    PrimerSpecificity,
    TargetContext,
    Workflow,
)
from smfprimer.output import outcome_rows
from smfprimer.specificity import assess_specificity, expand_degenerate


def _primer(sequence: str, side: str, start: int, end: int) -> CandidatePrimer:
    return CandidatePrimer(
        sequence=sequence,
        unconverted_sequence=sequence.replace("Y", "C").replace("R", "G"),
        side=side,
        start=start,
        end=end,
        tm=50,
        gc_fraction=0.5,
        degeneracies=sequence.count("Y") + sequence.count("R"),
        score=0,
    )


def _outcome() -> DesignOutcome:
    parameters = DesignParameters(
        min_length=2,
        optimum_length=2,
        max_length=2,
        min_tm=0,
        optimum_tm=10,
        max_tm=100,
        min_gc=0,
        max_gc=1,
        search_window=2,
        min_amplicon_size=10,
        max_amplicon_size=30,
    )
    target = DesignTarget("target", "intended", "A" * 40, 15, 16)
    pair = PrimerPair(
        _primer("AY", "forward", 5, 7),
        _primer("TR", "reverse", 23, 25),
        5,
        25,
        ConvertedStrand.TOP,
        0,
    )
    return DesignOutcome(
        target,
        Workflow.CONVERSION,
        TargetContext.BOTH,
        ConvertedStrand.TOP,
        parameters,
        (pair,),
    )


def _fake_bowtie(path: Path) -> Path:
    script = path / "fake-bowtie"
    script.write_text(
        "#!/usr/bin/env python3\n"
        "import sys\n"
        "from pathlib import Path\n"
        "records = []\n"
        "name = None\n"
        "for line in Path(sys.argv[-2]).read_text().splitlines():\n"
        "    if line.startswith('>'):\n"
        "        name = line[1:]\n"
        "    else:\n"
        "        records.append((name, line))\n"
        "lines = []\n"
        "for name, sequence in records:\n"
        "    reverse = '_reverse_' in name\n"
        "    flag = 16 if reverse else 0\n"
        "    position = 30 if reverse else 10\n"
        "    fields = [name, str(flag), 'chr1', str(position), '255',\n"
        "              f'{len(sequence)}M', '*', '0', '0', sequence, '*', 'NM:i:0']\n"
        "    lines.append('\\t'.join(fields))\n"
        "Path(sys.argv[-1]).write_text('\\n'.join(lines) + '\\n')\n"
    )
    script.chmod(0o755)
    return script


def test_expand_degenerate_iupac_sequence() -> None:
    assert expand_degenerate("ARY") == ("AAC", "AAT", "AGC", "AGT")


def test_bowtie_hits_are_deduplicated_and_paired(tmp_path: Path) -> None:
    outcomes = assess_specificity(
        [_outcome()],
        "unused-index",
        bowtie=str(_fake_bowtie(tmp_path)),
    )
    specificity = outcomes[0].pairs[0].specificity
    assert specificity is not None
    assert specificity.forward_hit_count == 1
    assert specificity.reverse_hit_count == 1
    assert specificity.intended_amplicon_count == 0
    assert specificity.off_target_amplicon_count == 1
    assert specificity.off_target_amplicons[0].reference_name == "chr1"
    assert specificity.off_target_amplicons[0].start == 9
    assert specificity.off_target_amplicons[0].end == 31

    row = outcome_rows(outcomes)[0]
    assert row["specificity_status"] == "off_targets"
    assert row["specificity_off_target_amplicons"] == 1
    assert json.loads(str(row["specificity_off_target_loci"]))[0]["length"] == 22


def test_status_requires_the_intended_product_to_call_a_pair_specific() -> None:
    no_products = PrimerSpecificity("index", 2, 1, 0, 0, ())
    intended_only = PrimerSpecificity("index", 2, 1, 1, 1, ())
    assert no_products.status == "no_off_targets"
    assert intended_only.status == "specific"
