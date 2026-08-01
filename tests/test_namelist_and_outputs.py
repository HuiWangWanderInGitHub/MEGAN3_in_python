#!/usr/bin/env python3
"""Tests for namelist-driven runs and unit-bearing output headers."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import tempfile
import unittest

import pandas as pd

from src.megan_model import load_inputs, load_namelist, simulate, write_outputs


PROJECT_ROOT = Path(__file__).resolve().parents[1]


class NamelistAndOutputTests(unittest.TestCase):
    """Verify that normal runs need only one configuration file."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.configuration = load_namelist(PROJECT_ROOT / "MEGAN_namelist.ini")
        cls.bundle = load_inputs(
            cls.configuration.paths,
            cls.configuration.controls,
        )
        cls.result = simulate(cls.bundle, cls.configuration.model)

    def test_namelist_contains_paths_controls_model_and_output_settings(self) -> None:
        """The namelist should fully define a normal model run."""

        configuration = self.configuration
        self.assertEqual(configuration.controls.humidity_mode, "RH")
        self.assertAlmostEqual(configuration.controls.latitude_deg, 45.4)
        self.assertEqual(configuration.model.number_of_layers, 5)
        self.assertTrue(configuration.output.create_diagnostic_figure)
        self.assertTrue(configuration.paths.meteorology.is_absolute())
        self.assertTrue(configuration.paths.output_directory.is_absolute())

    def test_written_csv_headers_include_units_and_namelist_is_copied(self) -> None:
        """Every CSV field, especially all species, must expose its unit."""

        with tempfile.TemporaryDirectory() as temporary_directory:
            output_directory = Path(temporary_directory)
            settings = replace(
                self.configuration.output,
                create_diagnostic_figure=False,
                show_diagnostic_figure=False,
            )
            paths = write_outputs(
                self.result,
                output_directory,
                output_options=settings,
                namelist_path=self.configuration.source_path,
            )

            all_columns = pd.read_csv(paths["all_species"], nrows=0).columns.tolist()
            isoprene_columns = pd.read_csv(paths["isoprene"], nrows=0).columns.tolist()

            self.assertEqual(all_columns[0], "day [day of year]")
            self.assertEqual(all_columns[1], "hour [local decimal hour]")
            self.assertEqual(len(all_columns[2:]), 19)
            self.assertTrue(
                all(column.endswith(" [nmol m-2 s-1]") for column in all_columns[2:])
            )
            self.assertIn("isoprene observed [mg m-2 h-1]", isoprene_columns)
            self.assertIn("air temperature [deg C]", isoprene_columns)
            self.assertIn("LAI [m2 m-2]", isoprene_columns)
            self.assertEqual(
                paths["namelist"].read_text(encoding="utf-8"),
                self.configuration.source_path.read_text(encoding="utf-8"),
            )


if __name__ == "__main__":
    unittest.main()
