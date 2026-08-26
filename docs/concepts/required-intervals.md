# Required intervals

A required interval is the span that every returned PCR product must contain.
Primer sites are placed outside the interval, and the complete amplicon must
satisfy `--product-size`.

## How each input mode defines it

| Mode | Required interval |
| --- | --- |
| `sequence` | `--target-start` through `--target-end` |
| `fasta` | Centered interval of `--target-width` bases |
| `annotated` | Every matching named feature |
| `bed` | Each BED row |
| `tss` | Requested upstream/downstream span around each TSS |

For GenBank and SnapGene input, a feature labeled `required_interval` is used
by default. Multiple such features produce independent design results. Other
labels can be selected with repeated `--required-feature` arguments or with
`--required-features-file`.

## Search window versus product size

`--search-window` controls how far smfprimer supplies sequence to Primer3 on
each side of the target. `--product-size` constrains the complete product.
Both must be large enough to fit the required interval and two primers.

If no pair is possible, the target is retained in output with
`status=no_candidates` rather than silently discarded.
