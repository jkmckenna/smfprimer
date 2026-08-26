# Design from GenBank or SnapGene

Annotated mode is useful for plasmids and constructs because existing labels
can define required intervals and remain visible in output.

## Default `required_interval`

The example GenBank file contains a feature labeled `required_interval`:

```console
smfprimer design annotated docs/examples/annotated_construct.gb \
  --workflow conversion \
  --product-size 50:100 \
  --min-tm 30 \
  --optimum-tm 50 \
  --max-tm 80 \
  --output annotated_primers.tsv
```

The companion `annotated_primers.gb` retains the original `example_promoter`
and `required_interval` features and adds ranked amplicon and primer features.

## Select custom feature names

Repeat the option to design around several labels:

```console
smfprimer design annotated construct.dna \
  --required-feature promoter_A \
  --required-feature enhancer_B \
  --product-size 300:700 \
  --output construct_primers.tsv
```

Or provide one label per line:

```console
smfprimer design annotated construct.gb \
  --required-features-file docs/examples/required_features.txt \
  --output construct_primers.tsv
```

Each matching feature becomes a separate target. SnapGene `.dna` files are
accepted as input; the companion annotation file is GenBank because Biopython
does not provide a maintained SnapGene writer.

!!! warning
    Required features that cross the origin of a circular construct are not
    currently supported. Represent the desired target as one contiguous span.
