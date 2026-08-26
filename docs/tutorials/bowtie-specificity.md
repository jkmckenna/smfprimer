# Assess specificity with Bowtie 1

Specificity assessment aligns concrete expansions of each ordered primer to a
Bowtie 1 index and identifies inward-facing products within the configured
amplicon-size range.

## Build an index

```console
bowtie-build GRCh38.fa GRCh38
```

The second argument is the index prefix passed to smfprimer.

## Add specificity to design

```console
smfprimer design fasta loci.fa \
  --product-size 250:500 \
  --bowtie-index GRCh38 \
  --specificity-mismatches 2 \
  --specificity-max-expansions 256 \
  --output primers.tsv
```

The tool invokes Bowtie 1 only when `--bowtie-index` is supplied. Use
`--bowtie-executable /path/to/bowtie` when it is not on `PATH`.

## Interpret the columns

- `specificity_forward_hits` and `specificity_reverse_hits` count individual
  primer alignments.
- `specificity_off_target_amplicons` counts compatible background products,
  which is more informative than individual hits alone.
- `specificity_off_target_loci` records coordinates and mismatch counts.
- `specificity_status=no_off_targets` means no compatible background product
  was formed under the configured mismatch and size rules.

Degenerate primers are expanded before alignment. A primer with many IUPAC
positions can exceed `--specificity-max-expansions`; increase the cap carefully
because alignment work grows with the number of expansions.

!!! important
    Specificity is defined by the supplied index. Use an index representing the
    actual experimental background, including an appropriate transformed
    reference when assessing converted DNA.
