# smfprimer

`smfprimer` designs conversion-aware PCR primer pairs for single-molecule
footprinting experiments. Primer3 selects and ranks the pairs; smfprimer adds
conversion chemistry, required-interval handling, annotated sequence output,
and optional Bowtie 1 specificity assessment.

## Choose a workflow

| Starting material | Start here |
| --- | --- |
| Amplicon or locus FASTA | [Amplicon FASTA tutorial](tutorials/fasta-input.md) |
| Annotated GenBank or SnapGene file | [Annotated-input tutorial](tutorials/annotated-input.md) |
| Genome FASTA plus BED/GTF | [Genome, BED, and GTF tutorial](tutorials/genome-bed-gtf.md) |
| Existing forward/reverse primers | [Primer evaluation tutorial](tutorials/evaluate-primers.md) |
| Background-genome specificity | [Bowtie specificity tutorial](tutorials/bowtie-specificity.md) |

## What a design produces

A successful file-based run creates two complementary artifacts:

- A TSV or JSON table containing ranked primer sequences, Primer3 metrics,
  coordinates, conversion ambiguity, and optional specificity metrics.
- A GenBank file containing the input annotations, required interval, ranked
  amplicon spans, and strand-aware primer annotations.

Coordinates in CLI arguments, BED input, TSV/JSON output, and the Python API
are zero-based and half-open. GTF and GenBank coordinates are converted at the
input/output boundary.

## Minimal example

```console
smfprimer design fasta docs/examples/amplicon.fa \
  --target-width 10 \
  --product-size 50:100 \
  --workflow conversion \
  --output primers.tsv
```

This writes `primers.tsv` and `primers.gb`.

!!! note
    Primer design is an assay-development aid. Review ambiguity, specificity,
    product length, and the annotated amplicon before ordering oligos.
