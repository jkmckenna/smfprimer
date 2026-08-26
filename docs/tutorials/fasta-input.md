# Design from an amplicon FASTA

Use FASTA mode when each record contains a locus or expected amplicon and the
required interval is centered within that record.

## Run the example

From the repository root:

```console
smfprimer design fasta docs/examples/amplicon.fa \
  --target-width 10 \
  --workflow conversion \
  --context both \
  --converted-strand top \
  --product-size 50:100 \
  --min-tm 30 \
  --optimum-tm 50 \
  --max-tm 80 \
  --output primers.tsv
```

The relaxed Tm range makes the short repetitive documentation fixture useful;
use assay-appropriate values for real designs.

## Inspect the result

Open `primers.tsv` and begin with:

- `status`: `ok` means Primer3 returned a pair.
- `rank`: lower ranks have lower Primer3 pair penalties.
- `forward_order_sequence` and `reverse_order_sequence`: final 5′→3′ oligos.
- `forward_degeneracies` and `reverse_degeneracies`: IUPAC ambiguity counts.
- `amplicon_start`, `amplicon_end`, and `amplicon_length`: predicted product.

Open `primers.gb` in a GenBank viewer to see the required interval, amplicon,
and primer directions together.

## Multiple records

Every FASTA record produces its own outcome. Failed records remain in the
table, making a batch run auditable.

Add `--gtf annotations.gtf.gz` when FASTA record names and coordinates match
the supplied annotation. Overlapping GTF features are included in the GenBank
companion.
