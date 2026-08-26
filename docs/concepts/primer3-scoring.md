# Primer3 scoring

Primer3 performs candidate selection and pair ranking. smfprimer maps its
design parameters directly to Primer3 size, Tm, GC, and product-size settings.

## Reading scores

`pair_score` is Primer3's pair penalty. Lower values rank ahead of higher
values. The individual forward and reverse primer penalties are retained in
the Python data model and in GenBank primer annotations.

The penalty is an objective function, not a probability or percentage. It
reflects departures from preferred values and Primer3's pair constraints.
Scores are most useful for comparing pairs from the same run with the same
parameters.

## Existing primers

`evaluate_primer_pairs` invokes Primer3's `check_primers` task. Supplied
primers are retained even when they violate a configured threshold, and the
corresponding Primer3 problems are returned as warnings.

## Conversion-aware limitation

Primer3 scores the concrete expected converted allele rather than every member
of a degenerate oligo mixture. Always inspect ambiguity counts. Highly
degenerate primers can also exceed the Bowtie expansion cap used during
specificity assessment.

`--max-degeneracies` (default: 3) discards candidate primers whose ambiguous
`Y`/`R` base count exceeds the threshold before pairs are ranked, since
Primer3 itself has no notion of ambiguity. The CLI accepts any non-negative
integer; to disable the filter entirely, use the Python API with
`DesignParameters(max_degeneracies=None)`.
