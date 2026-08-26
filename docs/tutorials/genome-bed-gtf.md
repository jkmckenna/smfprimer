# Design from genome, BED, and GTF files

BED mode treats every BED row as a required interval. GTF annotation is
optional but useful for preserving nearby genes and regulatory features in the
GenBank output.

## BED example

```console
smfprimer design bed \
  --genome docs/examples/genome.fa \
  --bed docs/examples/targets.bed \
  --gtf docs/examples/annotations.gtf \
  --search-window 40 \
  --product-size 50:100 \
  --min-tm 30 \
  --optimum-tm 50 \
  --max-tm 80 \
  --output genome_primers.tsv
```

BED and output coordinates are zero-based and half-open. GTF coordinates are
one-based and closed and are converted internally.

The genome FASTA must be uncompressed and have consistent sequence-line width
within each record. smfprimer builds an in-memory access index; it does not
write an index beside the reference.

## TSS example

TSS mode accepts a file containing one `gene_id` or `gene_name` per line:

```console
smfprimer design tss \
  --genome docs/examples/genome.fa \
  --gtf docs/examples/annotations.gtf \
  --genes docs/examples/genes.txt \
  --tss-upstream 5 \
  --tss-downstream 5 \
  --tss-policy all \
  --search-window 40 \
  --product-size 50:100 \
  --min-tm 30 \
  --optimum-tm 50 \
  --max-tm 80 \
  --output tss_primers.tsv
```

Upstream and downstream follow transcriptional orientation. Use
`--tss-policy longest` to choose the longest transcript per requested gene.
