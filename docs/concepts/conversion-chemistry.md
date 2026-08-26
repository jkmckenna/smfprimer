# Conversion chemistry

smfprimer always accepts the ordinary, unconverted top-strand reference. The
`--converted-strand` option identifies which reference strand undergoes the
experimental chemistry.

## Deaminase workflow

With `--workflow deaminase`, every cytosine on the converted strand is treated
as a possible C/T position:

- A forward ordered primer reports `Y` for C/T.
- A reverse ordered primer reports the complementary `R` for A/G.

## Conversion workflow

With `--workflow conversion`, non-target cytosines become expected thymine.
Cytosines in the selected protected context remain ambiguous:

- `--context cpg` protects CpG cytosines.
- `--context gpc` protects GpC cytosines.
- `--context both` protects both groups.

## How Primer3 sees the sequence

Primer3 requires a concrete design sequence. smfprimer uses the expected
converted member of each mixture (`Y` as `T`, and its ordered-primer complement
`R` as `A`) for Primer3 selection and thermodynamic metrics. It then restores
the conversion-aware IUPAC bases at Primer3's selected coordinates.

Consequently:

- `forward_order_sequence` and `reverse_order_sequence` are the sequences to
  review for ordering.
- Tm, GC, and Primer3 penalties describe the concrete converted allele.
- `forward_degeneracies` and `reverse_degeneracies` show how many ambiguous
  positions each ordered mixture contains.

Do not interpret a low Primer3 penalty as a specificity guarantee. Specificity
is assessed independently against the chosen Bowtie index.
