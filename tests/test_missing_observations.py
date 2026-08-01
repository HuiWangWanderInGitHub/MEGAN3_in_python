#!/usr/bin/env python3
"""Tests for optional and partially missing isoprene observations."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

from src.megan_model import load_inputs, load_namelist, simulate


PROJECT_ROOT = Path(__file__).resolve().parents[1]


class MissingIsopreneObservationTests(unittest.TestCase):
    """Observation gaps must not affect the modeled emissions."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.configuration = load_namelist(PROJECT_ROOT / "MEGAN_namelist.ini")
        cls.original_frame = pd.read_csv(
            cls.configuration.paths.meteorology,
            encoding="utf-8-sig",
        )
        cls.reference = simulate(
            load_inputs(cls.configuration.paths, cls.configuration.controls),
            cls.configuration.model,
        )

    def _simulate_with_frame(self, frame: pd.DataFrame):
        """Write a temporary meteorology file and run the normal model path."""

        with tempfile.TemporaryDirectory() as temporary_directory:
            meteorology_path = Path(temporary_directory) / "meteorology.csv"
            frame.to_csv(meteorology_path, index=False, na_rep="")
            paths = replace(
                self.configuration.paths,
                meteorology=meteorology_path,
                output_directory=Path(temporary_directory) / "outputs",
            )
            bundle = load_inputs(paths, self.configuration.controls)
            result = simulate(bundle, self.configuration.model)
            return bundle, result

    def test_blank_observation_cells_are_skipped_only_in_diagnostics(self) -> None:
        """Blank Isop_obs cells become NaN without changing model emissions."""

        frame = self.original_frame.copy()
        observation_column = "Isop(mg/m2/h)"
        daytime_finite = (
            frame["Hour"].between(9.0, 17.0, inclusive="both")
            & frame[observation_column].notna()
        )
        rows_to_blank = frame.index[daytime_finite][:5]
        self.assertEqual(len(rows_to_blank), 5)
        frame.loc[rows_to_blank, observation_column] = np.nan

        bundle, result = self._simulate_with_frame(frame)

        self.assertTrue(
            bundle.meteorology.loc[rows_to_blank, "isoprene_observed_mg_m2_h"]
            .isna()
            .all()
        )
        np.testing.assert_allclose(
            result.all_species.iloc[:, 2:].to_numpy(dtype=float),
            self.reference.all_species.iloc[:, 2:].to_numpy(dtype=float),
            rtol=0.0,
            atol=0.0,
            equal_nan=True,
        )
        self.assertEqual(
            result.metrics["n_observed_available"],
            self.reference.metrics["n_observed_available"] - 5,
        )
        self.assertEqual(
            result.metrics["n"],
            self.reference.metrics["n"] - 5,
        )

    def test_observation_column_may_be_omitted_entirely(self) -> None:
        """A simulation can run when no isoprene observation column exists."""

        frame = self.original_frame.drop(columns=["Isop(mg/m2/h)"])
        bundle, result = self._simulate_with_frame(frame)

        self.assertTrue(
            bundle.meteorology["isoprene_observed_mg_m2_h"].isna().all()
        )
        self.assertEqual(result.metrics["n_observed_available"], 0)
        self.assertEqual(result.metrics["n"], 0)
        np.testing.assert_allclose(
            result.all_species.iloc[:, 2:].to_numpy(dtype=float),
            self.reference.all_species.iloc[:, 2:].to_numpy(dtype=float),
            rtol=0.0,
            atol=0.0,
            equal_nan=True,
        )


if __name__ == "__main__":
    unittest.main()
