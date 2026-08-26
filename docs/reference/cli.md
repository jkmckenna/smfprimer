# Command-line reference

The command hierarchy is:

```text
smfprimer design {sequence,fasta,annotated,bed,tss} [options]
```

Run `smfprimer design MODE --help` for the authoritative option list and
defaults installed with your version.

## Shared design options

| Option | Meaning |
| --- | --- |
| `--workflow {deaminase,conversion}` | Experimental conversion model |
| `--converted-strand {top,bottom}` | Reference strand subjected to conversion |
| `--context {cpg,gpc,both}` | Protected cytosine context |
| `--min-length`, `--optimum-length`, `--max-length` | Primer3 size constraints |
| `--min-tm`, `--optimum-tm`, `--max-tm` | Primer3 Tm constraints |
| `--min-gc`, `--max-gc` | GC fractions from 0 to 1 |
| `--search-window N` | Bases searched on each side of the target |
| `--product-size MIN:MAX` | Complete amplicon-size range |
| `--max-results N` | Maximum ranked pairs per target |
| `--max-degeneracies N` | Maximum degenerate (`Y`/`R`) bases tolerated per primer (default: 3) |
| `--format {tsv,json}` | Tabular output format |
| `--output PATH` | Write table rather than standard output |
| `--genbank-output PATH` | Override companion GenBank path |
| `--no-genbank` | Disable companion GenBank output |
| `--bowtie-index PREFIX` | Enable Bowtie 1 specificity |

## Input-specific options

### Sequence

```console
smfprimer design sequence SEQUENCE_OR_FASTA \
  --target-start START --target-end END
```

`--target-id` controls the result name.

### FASTA

```console
smfprimer design fasta LOCI.fa --target-width N [--gtf FEATURES.gtf.gz]
```

### Annotated

```console
smfprimer design annotated CONSTRUCT.gb \
  [--required-feature LABEL] \
  [--required-features-file NAMES.txt]
```

### BED

```console
smfprimer design bed --genome GENOME.fa --bed TARGETS.bed [--gtf FEATURES.gtf.gz]
```

### TSS

```console
smfprimer design tss \
  --genome GENOME.fa --gtf FEATURES.gtf --genes GENES.txt \
  [--tss-upstream N] [--tss-downstream N] \
  [--tss-policy {all,longest}]
```
