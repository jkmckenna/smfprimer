# Python API

## Design one interval

```python
from smfprimer import DesignParameters, Workflow, design_primers

pairs = design_primers(
    reference_sequence,
    target_start=250,
    target_end=450,
    workflow=Workflow.CONVERSION,
    parameters=DesignParameters(
        min_amplicon_size=300,
        max_amplicon_size=700,
    ),
)
```

Each `PrimerPair` contains forward and reverse `CandidatePrimer` records,
top-reference amplicon coordinates, converted-strand identity, Primer3 pair
penalty, and optional specificity.

## Design normalized targets

```python
from smfprimer import DesignParameters, TargetContext, design_targets
from smfprimer.targets import bed_targets

targets = bed_targets("GRCh38.fa", "targets.bed", flank=500)
outcomes = design_targets(
    targets,
    context=TargetContext.BOTH,
    parameters=DesignParameters(
        min_amplicon_size=300,
        max_amplicon_size=700,
    ),
)
```

`DesignOutcome` preserves targets with no returned pairs, along with a status
and explanation.

## Evaluate supplied pairs

```python
from smfprimer import PrimerSet, evaluate_primer_pairs

results = evaluate_primer_pairs(
    [
        PrimerSet("pair_1", "ACGT...", "TGCA..."),
        PrimerSet("pair_2", "GCTA...", "TAGC..."),
    ],
    reference_sequence,
    template_name="locus",
)
```

Supply `bowtie_index="GRCh38"` to attach specificity metrics.

## Serialize outcomes

```python
from smfprimer import format_genbank
from smfprimer.output import format_json, format_tsv

tsv_text = format_tsv(outcomes)
json_text = format_json(outcomes)
genbank_text = format_genbank(outcomes)
```
