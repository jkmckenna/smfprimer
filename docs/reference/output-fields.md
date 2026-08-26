# Output fields

TSV and JSON contain the same logical fields. Empty pair-specific fields on a
`no_candidates` row are expected.

## Target and run

| Field | Meaning |
| --- | --- |
| `target_id` | Stable target label |
| `source` | Input adapter, such as `fasta`, `bed`, or `snapgene` |
| `reference_name` | Source reference identifier |
| `feature_strand` | Strand of the source annotation when available |
| `target_start`, `target_end` | Required interval on the source reference |
| `status`, `message` | Design outcome and explanation |
| `rank` | Primer3 pair rank, starting at 1 |
| `workflow`, `context`, `converted_strand` | Chemistry settings |
| `design_engine` | `primer3` |
| `metadata`, `parameters` | JSON objects describing source and design settings |

## Amplicon

| Field | Meaning |
| --- | --- |
| `amplicon_start`, `amplicon_end`, `amplicon_length` | Complete product coordinates and size |
| `amplicon_reference_sequence` | Original top-reference product sequence |
| `amplicon_unconverted_sequence` | Product in converted-strand orientation before conversion |
| `amplicon_converted_sequence` | Predicted conversion-aware product |

## Primers

Forward and reverse fields use the same suffixes:

| Suffix | Meaning |
| --- | --- |
| `_start`, `_end` | Binding-site coordinates on the top reference |
| `_unconverted_sequence` | Binding sequence before conversion, in primer orientation |
| `_order_sequence` | Final 5′→3′ IUPAC sequence to review for ordering |
| `_tm` | Primer3 Tm for the concrete expected converted allele |
| `_gc_fraction` | Primer3 GC percentage divided by 100 |
| `_degeneracies` | Number of non-ACGT bases in the order sequence |

`pair_score` is Primer3's pair penalty; lower ranks are preferred within the
same run.

## Specificity

Specificity fields are populated only when a Bowtie index is supplied:

- `specificity_status`
- `specificity_index`
- `specificity_mismatches`
- `specificity_forward_hits`
- `specificity_reverse_hits`
- `specificity_intended_amplicons`
- `specificity_off_target_amplicons`
- `specificity_off_target_loci`

The GenBank companion stores the same design context as feature qualifiers and
retains source annotations when the input supplies them.
