# smfprimer

`smfprimer` is a minimal Python package and command-line tool for designing
primers for experiments analyzed with
[smftools](https://github.com/jkmckenna/smftools).

It supports deaminase and conversion chemistries, explicit converted-strand
selection, and CpG-, GpC-, or dual-context conversion designs.

Tutorials, concepts, and complete CLI/output documentation live in the
[documentation site source](docs/index.md). Build it locally with
`mkdocs serve`.

## Chemistry model

- **Deaminase:** every cytosine on the converted strand is uncertain. Forward
  primers use `Y` (C/T), while reverse primers use the complementary `R` (A/G).
- **Conversion:** non-target cytosines are represented as converted thymine.
  Cytosines in selected CpG/GpC contexts can be protected, so reported oligos
  use `Y`/`R` at those uncertain positions.
- **Stranding:** the input is always a top-strand reference. Selecting bottom
  conversion reverses the chemistry internally, then maps all coordinates back
  to the supplied reference.

Primer3 selects and ranks primer pairs using its thermodynamic model, size, Tm,
GC, pair-complementarity, and product-size constraints. Design runs against the
expected concrete converted allele (`Y` as `T`); smfprimer then reconstructs
the order sequence from the IUPAC conversion template at Primer3's selected
coordinates. Reported Tm, GC, primer scores, and pair scores are Primer3 values
for that concrete allele. Ambiguity counts remain explicit for review.

## Installation

```console
python -m pip install -e .
```

The installed command has five input modes:

```console
smfprimer design {sequence,fasta,annotated,bed,tss} --help
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

Add `--gtf annotations.gtf.gz` to retain every overlapping GTF feature in the
annotated GenBank output. GTF coordinates are clipped to each FASTA record and
remapped to record-local coordinates.

## Annotated GenBank or SnapGene input

GenBank (`.gb`, `.gbk`, `.gbff`, `.genbank`) and native SnapGene (`.dna`)
inputs retain their original feature labels, qualifiers, coordinates, strands,
and topology. By default, every feature labeled `required_interval` becomes an
independent required primer-design interval:

```console
smfprimer design annotated construct.dna \
  --product-size 300:700 \
  --output primers.tsv
```

Multiple matching features produce multiple design outcomes. To use other
feature labels or names, repeat `--required-feature` or provide one name per
line with `--required-features-file`:

```console
smfprimer design annotated construct.gb \
  --required-feature promoter_A \
  --required-feature enhancer_B \
  --output primers.tsv

smfprimer design annotated construct.dna \
  --required-features-file required_features.txt \
  --output primers.tsv
```

SnapGene files are accepted as input and the annotated companion output is
GenBank. Native `.dna` writing is intentionally not attempted because
Biopython does not provide a maintained SnapGene writer.

## Genome FASTA plus BED

Every BED interval is a required span. Forward and reverse primer sites must
flank the entire interval, and the resulting product must satisfy the requested
size range.

```console
smfprimer design bed \
  --genome hg38.fa \
  --bed targets.bed \
  --gtf genes.gtf.gz \
  --product-size 300:700 \
  --search-window 500
```

Genome FASTA access is indexed in memory without loading chromosome sequences
or writing an index beside the reference. The FASTA must be uncompressed and
use a consistent sequence-line width within each record. When `--gtf` is
provided, overlapping features are clipped and retained in the GenBank record.

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
TSS-mode GenBank output automatically includes overlapping GTF features, and
both plain `.gtf` and `.gtf.gz` inputs are supported.

## Output

TSV is written to standard output by default. Use `--output PATH` to write a
file or `--format json` for structured JSON. File-based output automatically
gets a companion `.gb` GenBank file containing each design template with the
required target, ranked amplicons, and strand-aware primer annotations. Use
`--genbank-output PATH` to choose its location or `--no-genbank` to disable it.
Failed targets remain in batch output with `status=no_candidates` and an
explanation.

## Evaluate existing primer pairs

`evaluate_primer_pairs` accepts one `PrimerSet` or any iterable of sets. Primer3
`check_primers` supplies each primer's thermodynamic metrics, constraint
problems, and Primer3 penalty, plus the pair penalty. smfprimer adds compatible
template products, ambiguity counts, and optional Bowtie 1 specificity metrics:

```python
from smfprimer import PrimerSet, Workflow, evaluate_primer_pairs

results = evaluate_primer_pairs(
    [
        PrimerSet(
            "JM351_JM349",
            forward="GAGAGATYTGGYAGYGGAGAG",
            reverse="CTTTTCTRTCACCAATCCTRTCCC",
        )
    ],
    template_sequence,
    template_name="JM135",
    workflow=Workflow.DEAMINASE,
    bowtie_index="GRCh38",
    mismatches=2,
)
```

The template must be the unconverted top-strand sequence. Set
`converted_strand="bottom"` when evaluating primers for the converted bottom
strand. Product-size limits and primer thresholds can be supplied with a
`DesignParameters` instance.

## Bowtie 1 specificity assessment

Pass a Bowtie 1 index prefix to align every concrete expansion of each ordered
primer and report inward-facing genomic products within `--product-size`:

```console
bowtie-build GRCh38.fa GRCh38
smfprimer design fasta loci.fa \
  --product-size 250:500 \
  --bowtie-index GRCh38 \
  --specificity-mismatches 2 \
  --output primers.tsv
```

Specificity columns report primer hit counts, intended products, off-target
product counts, and off-target loci. Bowtie 1 is optional and is only invoked
when `--bowtie-index` is supplied. The index defines the exact background being
tested; use an appropriately transformed index when assessing converted-DNA
backgrounds.

Each successful pair includes:

- Top-reference target, primer-site, and amplicon coordinates
- The original top-reference amplicon
- The unconverted amplicon in converted-strand orientation
- The predicted converted amplicon
- Forward and reverse sequences before conversion, in primer orientation
- Final 5′→3′ `forward_order_sequence` and `reverse_order_sequence`
- Primer3 Tm, GC fraction, primer/pair penalties, ambiguity count, parameters,
  design-engine provenance, and source metadata

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
