#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Run the site-scale MEGAN model from one external namelist file.

Normal use
----------
Edit ``MEGAN_namelist.ini`` and run::

    python main_program.py

An alternative namelist can be selected without changing this program::

    python main_program.py --namelist experiments/site_b.ini

All input paths, model switches, numerical settings, diagnostics, output names,
and output formatting are defined in the namelist.  This file intentionally
contains no site-specific model parameters.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import time

from src.megan_model import (
    load_inputs,
    load_namelist,
    save_diagnostic_figure,
    simulate,
    write_outputs,
)


def build_parser() -> argparse.ArgumentParser:
    """Create a minimal parser containing only the namelist selector."""

    project_directory = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description=(
            "Run the site-scale MEGANv3.1 Python model using settings from "
            "a namelist-style INI file."
        )
    )
    parser.add_argument(
        "--namelist",
        type=Path,
        default=project_directory / "MEGAN_namelist.ini",
        help=(
            "Path to the run namelist. Relative paths inside the namelist are "
            "resolved relative to the namelist file itself."
        ),
    )
    return parser


def main() -> int:
    """Load the namelist, run the model, and write configured outputs."""

    arguments = build_parser().parse_args()
    configuration = load_namelist(arguments.namelist)

    print("MEGANv3.1 site model")
    print(f"  namelist: {configuration.source_path}")
    print(f"  meteorology: {configuration.paths.meteorology}")
    print(f"  output directory: {configuration.paths.output_directory}")

    start_time = time.perf_counter()
    inputs = load_inputs(configuration.paths, configuration.controls)
    print(
        "Inputs loaded: "
        f"{len(inputs.meteorology)} records, "
        f"{len(inputs.species.names)} species classes, "
        f"{sum(inputs.pfts.fractions_percent > 0)} active PFT(s)."
    )

    result = simulate(inputs, configuration.model)
    output_paths = write_outputs(
        result,
        configuration.paths.output_directory,
        output_options=configuration.output,
        namelist_path=configuration.source_path,
    )

    if configuration.output.create_diagnostic_figure:
        figure_path = (
            configuration.paths.output_directory
            / configuration.output.diagnostic_figure_filename
        )
        save_diagnostic_figure(
            result,
            figure_path,
            daytime_start_hour=configuration.model.daytime_start_hour,
            daytime_end_hour=configuration.model.daytime_end_hour,
            show=configuration.output.show_diagnostic_figure,
        )
        output_paths["figure"] = figure_path

    elapsed = time.perf_counter() - start_time
    print(f"Simulation finished in {elapsed:.2f} s.")
    for label, path in output_paths.items():
        print(f"  {label}: {path}")

    metrics = result.metrics
    print(
        "Daytime isoprene observations: "
        f"{metrics['n_observed_available']}/{metrics['n_daytime_records']} available; "
        f"{metrics['n']} paired finite records used."
    )
    if metrics["n"] and metrics["slope"] is not None:
        print(
            "Daytime isoprene diagnostics: "
            f"n={metrics['n']}, "
            f"slope={metrics['slope']:.3f}, "
            f"intercept={metrics['intercept']:.3f}, "
            f"R2={metrics['r_squared']:.3f}, "
            f"RMSE={metrics['rmse_mg_m2_h']:.3f} mg m-2 h-1."
        )
    elif metrics["n"]:
        print(
            "Daytime isoprene diagnostics: "
            f"n={metrics['n']}, "
            f"RMSE={metrics['rmse_mg_m2_h']:.3f} mg m-2 h-1; "
            "regression unavailable because the sample is too small or constant."
        )
    else:
        print("No paired daytime isoprene observations were available.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
