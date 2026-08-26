# Evaluate existing primer pairs

Use the Python API when primers already exist and you want Primer3 metrics,
template-product validation, and optional background specificity.

```python
from Bio import SeqIO

from smfprimer import (
    DesignParameters,
    PrimerSet,
    Workflow,
    evaluate_primer_pairs,
)

record = SeqIO.read("docs/examples/amplicon.fa", "fasta")

result = evaluate_primer_pairs(
    PrimerSet(
        name="example_pair",
        forward="GTGGGYGATTTTAGTAATTATT",
        reverse="CCATATCAARCTACRATATCTC",
    ),
    str(record.seq),
    template_name=record.id,
    workflow=Workflow.CONVERSION,
    parameters=DesignParameters(
        min_tm=30,
        optimum_tm=50,
        max_tm=80,
        min_amplicon_size=50,
        max_amplicon_size=100,
    ),
)[0]

print(result.status)
print(result.forward.tm, result.reverse.tm)
print(result.pair_score)
for product in result.valid_template_amplicons:
    print(product.start, product.end, product.length)
```

`PrimerSet` also accepts IUPAC bases such as `Y` and `R`. Primer3 metrics are
calculated for the concrete expected converted member of the mixture, while
template matching remains IUPAC-aware.

## Status values

- `pass`: primer constraints pass and a valid template product exists.
- `constraint_warning`: Primer3 reports a supplied-primer problem.
- `no_template_product`: no inward-facing product satisfies the size range.
- `off_targets`: Bowtie found one or more background products.

Inspect `forward.warnings`, `reverse.warnings`, and pair-level `warnings` rather
than relying on status alone.
