# smfprimer

`smfprimer` is a minimal Python package and command-line tool for designing
primers for experiments analyzed with
[smftools](https://github.com/jkmckenna/smftools).

It supports deaminase and conversion chemistries, explicit converted-strand
selection, and CpG-, GpC-, or dual-context conversion designs.

## Chemistry model

- **Deaminase:** every cytosine on the converted strand is uncertain. Forward
  primers use `Y` (C/T), while reverse primers use the complementary `R` (A/G).
- **Conversion:** non-target cytosines are represented as converted thymine.
  Cytosines in selected CpG/GpC contexts can be protected, so candidates avoid
  those positions when possible and use `Y`/`R` only as a fallback.
- **Stranding:** the input is always a top-strand reference. Selecting bottom
  conversion reverses the chemistry internally, then maps all coordinates back
  to the supplied reference.

Candidates without ambiguous bases rank ahead of candidates containing them.
The current version uses a dependency-free Wallace-rule Tm estimate for
screening; final thermodynamic and specificity checks are still required before
ordering primers.

## Installation

```console
python -m pip install -e .
```

The installed command has four input modes:

```console
smfprimer design {sequence,fasta,bed,tss} --help
```

Coordinates are zero-based and half-open in BED, CLI arguments, and output.
GTF's one-based closed coordinates are converted internally.

## One sequence

Provide DNA directly or a single-record FASTA:

```console
smfprimer design sequence AACCGG... \
  --target-start 250 \
  --target-end 450 \
  --product-size 300:700 \
  --workflow conversion \
  --context both \
  --converted-strand top
```

## Multi-FASTA

One target is created at the center of every record. `--target-width` controls
the required centered span; the complete product must satisfy
`--product-size`.

```console
smfprimer design fasta loci.fa \
  --target-width 25 \
  --product-size 250:500 \
  --output primers.tsv
```

## Genome FASTA plus BED

Every BED interval is a required span. Forward and reverse primer sites must
flank the entire interval, and the resulting product must satisfy the requested
size range.

```console
smfprimer design bed \
  --genome hg38.fa \
  --bed targets.bed \
  --product-size 300:700 \
  --search-window 500
```

Genome FASTA access is indexed in memory without loading chromosome sequences
or writing an index beside the reference. The FASTA must be uncompressed and
use a consistent sequence-line width within each record.

## TSS-centered targets

Gene lists contain one `gene_id` or `gene_name` per line. By default, one
target is produced for every unique transcript TSS. `--tss-policy longest`
instead selects the longest transcript for each requested gene.

```console
smfprimer design tss \
  --genome hg38.fa \
  --gtf genes.gtf \
  --genes genes.txt \
  --tss-upstream 200 \
  --tss-downstream 200 \
  --tss-policy all \
  --product-size 450:800
```

Upstream and downstream distances follow transcriptional orientation.

## Output

TSV is written to standard output by default. Use `--output PATH` to write a
file or `--format json` for structured JSON. Failed targets remain in batch
output with `status=no_candidates` and an explanation.

Each successful pair includes:

- Top-reference target, primer-site, and amplicon coordinates
- The original top-reference amplicon
- The unconverted amplicon in converted-strand orientation
- The predicted converted amplicon
- Forward and reverse sequences before conversion, in primer orientation
- Final 5′→3′ `forward_order_sequence` and `reverse_order_sequence`
- Tm, GC fraction, ambiguity count, pair score, parameters, and source metadata

This makes the reference sequence auditable while keeping the sequences that
should be ordered explicit.

## Python API

```python
from smfprimer import DesignParameters, TargetContext, design_targets
from smfprimer.targets import bed_targets

parameters = DesignParameters(
    min_amplicon_size=300,
    max_amplicon_size=700,
)
targets = bed_targets("hg38.fa", "targets.bed", flank=500)
pairs = design_targets(
    targets,
    context=TargetContext.BOTH,
    parameters=parameters,
)
```
