"""Command-line interface for smfprimer."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import NamedTuple

from .design import design_targets
from .models import ConvertedStrand, DesignParameters, DesignTarget, TargetContext, Workflow
from .output import format_json, format_tsv
from .targets import bed_targets, fasta_targets, read_fasta, sequence_target, tss_targets


class _ProductSize(NamedTuple):
    minimum: int
    maximum: int

    def __str__(self) -> str:
        return f"{self.minimum}:{self.maximum}"


def _product_size(value: str) -> _ProductSize:
    try:
        minimum_text, maximum_text = value.split(":", 1)
        minimum, maximum = int(minimum_text), int(maximum_text)
    except ValueError as error:
        raise argparse.ArgumentTypeError("product size must have the form MIN:MAX") from error
    if minimum < 1 or maximum < minimum:
        raise argparse.ArgumentTypeError("product size must satisfy 1 <= MIN <= MAX")
    return _ProductSize(minimum, maximum)


def _common_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    defaults = DesignParameters()
    parser.add_argument(
        "--workflow",
        choices=list(Workflow),
        default=Workflow.CONVERSION,
        help="experimental primer chemistry",
    )
    parser.add_argument(
        "--converted-strand",
        "--strand",
        choices=list(ConvertedStrand),
        default=ConvertedStrand.TOP,
        help="reference strand subjected to conversion",
    )
    parser.add_argument(
        "--context",
        choices=list(TargetContext),
        default=TargetContext.BOTH,
        help="protected conversion context",
    )
    parser.add_argument(
        "--min-length", type=int, default=defaults.min_length, help="minimum primer length"
    )
    parser.add_argument(
        "--optimum-length",
        type=int,
        default=defaults.optimum_length,
        help="preferred primer length",
    )
    parser.add_argument(
        "--max-length", type=int, default=defaults.max_length, help="maximum primer length"
    )
    parser.add_argument(
        "--min-tm", type=float, default=defaults.min_tm, help="minimum estimated primer Tm"
    )
    parser.add_argument(
        "--optimum-tm",
        type=float,
        default=defaults.optimum_tm,
        help="preferred estimated primer Tm",
    )
    parser.add_argument(
        "--max-tm", type=float, default=defaults.max_tm, help="maximum estimated primer Tm"
    )
    parser.add_argument(
        "--min-gc", type=float, default=defaults.min_gc, help="minimum primer GC fraction"
    )
    parser.add_argument(
        "--max-gc", type=float, default=defaults.max_gc, help="maximum primer GC fraction"
    )
    parser.add_argument(
        "--search-window",
        type=int,
        default=defaults.search_window,
        help="maximum primer search distance from each target edge",
    )
    parser.add_argument(
        "--max-results",
        type=int,
        default=defaults.max_results,
        help="maximum pairs returned per target",
    )
    parser.add_argument(
        "--product-size",
        type=_product_size,
        default=_ProductSize(100, 500),
        metavar="MIN:MAX",
        help="allowed complete amplicon length",
    )
    parser.add_argument("--format", choices=["tsv", "json"], default="tsv")
    parser.add_argument("--output", type=Path, help="output path; defaults to standard output")
    return parser


def _parser() -> argparse.ArgumentParser:
    formatter = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(
        prog="smfprimer",
        description=__doc__,
        formatter_class=formatter,
    )
    parser.add_argument("--version", action="version", version="%(prog)s 0.1.0")
    commands = parser.add_subparsers(dest="command", required=True)
    design = commands.add_parser("design", help="design primer pairs")
    modes = design.add_subparsers(dest="mode", required=True)
    common = _common_parser()

    sequence = modes.add_parser(
        "sequence",
        parents=[common],
        formatter_class=formatter,
        help="design around an interval in one sequence",
    )
    sequence.add_argument("reference", help="DNA sequence or single-record FASTA path")
    sequence.add_argument("--target-start", type=int, required=True, help="0-based target start")
    sequence.add_argument("--target-end", type=int, required=True, help="exclusive target end")
    sequence.add_argument("--target-id", default="target")

    fasta = modes.add_parser(
        "fasta",
        parents=[common],
        formatter_class=formatter,
        help="design one centered amplicon per FASTA record",
    )
    fasta.add_argument("fasta", type=Path, help="multi-FASTA design templates")
    fasta.add_argument(
        "--target-width",
        type=int,
        default=1,
        help="width of the required interval centered in each record",
    )

    bed = modes.add_parser(
        "bed",
        parents=[common],
        formatter_class=formatter,
        help="design around required BED intervals",
    )
    bed.add_argument("--genome", type=Path, required=True)
    bed.add_argument("--bed", type=Path, required=True)
    bed.add_argument(
        "--flank",
        type=int,
        help="reference bases loaded on each side; defaults to search-window",
    )

    tss = modes.add_parser(
        "tss",
        parents=[common],
        formatter_class=formatter,
        help="design TSS-centered amplicons from a GTF and gene list",
    )
    tss.add_argument("--genome", type=Path, required=True)
    tss.add_argument("--gtf", type=Path, required=True)
    tss.add_argument("--genes", type=Path, required=True)
    tss.add_argument("--tss-upstream", type=int, default=100)
    tss.add_argument("--tss-downstream", type=int, default=100)
    tss.add_argument("--tss-policy", choices=["all", "longest"], default="all")
    tss.add_argument(
        "--flank",
        type=int,
        help="reference bases loaded beyond the TSS interval; defaults to search-window",
    )
    return parser


def _single_sequence(value: str) -> str:
    path = Path(value)
    if not path.is_file():
        return value
    records = read_fasta(path)
    if len(records) != 1:
        raise ValueError("sequence mode requires a single-record FASTA")
    return next(iter(records.values()))


def _targets(args: argparse.Namespace) -> list[DesignTarget]:
    if args.mode == "sequence":
        return [
            sequence_target(
                _single_sequence(args.reference),
                args.target_start,
                args.target_end,
                target_id=args.target_id,
            )
        ]
    if args.mode == "fasta":
        return fasta_targets(args.fasta, target_width=args.target_width)
    flank = args.search_window if args.flank is None else args.flank
    if args.mode == "bed":
        return bed_targets(args.genome, args.bed, flank=flank)
    return tss_targets(
        args.genome,
        args.gtf,
        args.genes,
        upstream=args.tss_upstream,
        downstream=args.tss_downstream,
        flank=flank,
        policy=args.tss_policy,
    )


def main(argv: list[str] | None = None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)
    minimum_product, maximum_product = args.product_size
    try:
        parameters = DesignParameters(
            min_length=args.min_length,
            optimum_length=args.optimum_length,
            max_length=args.max_length,
            min_tm=args.min_tm,
            optimum_tm=args.optimum_tm,
            max_tm=args.max_tm,
            min_gc=args.min_gc,
            max_gc=args.max_gc,
            search_window=args.search_window,
            max_results=args.max_results,
            min_amplicon_size=minimum_product,
            max_amplicon_size=maximum_product,
        )
        outcomes = design_targets(
            _targets(args),
            workflow=args.workflow,
            converted_strand=args.converted_strand,
            context=args.context,
            parameters=parameters,
        )
    except (OSError, UnicodeError, ValueError) as error:
        parser.error(str(error))

    content = format_tsv(outcomes) if args.format == "tsv" else format_json(outcomes)
    if args.output:
        args.output.write_text(content)
    else:
        sys.stdout.write(content)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
