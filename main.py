"""Command-line interface for the HarmoClimate training pipeline."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parent / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from harmoclimate.display import DISPLAY_VARIABLE_CHOICES  # noqa: E402
from harmoclimate.pipeline import (  # noqa: E402
    clean_pipeline,
    display_pipeline,
    generate_pipeline,
    regenerate_pipeline,
    run_pipeline,
    template_pipeline,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run HarmoClimate model generation workflows.")
    subparsers = parser.add_subparsers(dest="command")

    generate_parser = subparsers.add_parser(
        "generate",
        help="Stream data for a station code and generate model artefacts.",
    )
    generate_parser.add_argument(
        "station_code",
        help="Météo-France station code (NUM_POSTE) used to filter the source dataset.",
    )

    regenerate_parser = subparsers.add_parser(
        "regenerate",
        help="Regenerate an existing model using cached data or by refetching the archives.",
    )
    regenerate_parser.add_argument(
        "model_json",
        help="Path or file name of a previously generated model JSON (e.g. fr_bourges_temperature.json).",
    )

    display_parser = subparsers.add_parser(
        "display",
        help="Render a plot for a target model JSON and save it under generated/media.",
    )
    display_parser.add_argument(
        "model_json",
        help="Path or file name of a generated model JSON (e.g. fr_bourges_temperature.json).",
    )
    display_parser.add_argument(
        "--mode",
        choices=("annual", "intraday", "stochastic_annual", "stochastic_intraday"),
        default="annual",
        help="Choose 'annual' or 'stochastic_annual' for yearly plots, 'intraday' or 'stochastic_intraday' for daily.",
    )
    display_parser.add_argument(
        "--day",
        type=int,
        help="Solar day number (1-based) required when --mode=intraday.",
    )
    display_parser.add_argument(
        "--variables",
        nargs="+",
        choices=DISPLAY_VARIABLE_CHOICES,
        metavar="VAR",
        help=(
            "Ordered list of panel codes to render (choose from T, RH, TD, Q, E, P). "
            "Defaults to T Q P when omitted."
        ),
    )

    template_parser = subparsers.add_parser(
        "template",
        help="Generate an embedded template for an existing station model.",
    )
    template_parser.add_argument(
        "model_name",
        help="Model basename or JSON file (e.g. fr_bourges or fr_bourges_temperature.json).",
    )
    template_parser.add_argument(
        "language",
        help="Target template language (currently only 'cpp').",
    )

    subparsers.add_parser(
        "clean",
        help="Remove cached Parquet datasets under generated/data.",
    )

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "generate":
        generate_pipeline(args.station_code)
        return 0
    if args.command == "regenerate":
        regenerate_pipeline(args.model_json)
        return 0
    if args.command == "display":
        display_pipeline(
            args.model_json,
            mode=args.mode,
            day=args.day,
            variables=args.variables,
        )
        return 0
    if args.command == "template":
        template_pipeline(args.model_name, args.language)
        return 0
    if args.command == "clean":
        clean_pipeline()
        return 0

    # Backwards-compatible behaviour: default to the configured station.
    run_pipeline()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
