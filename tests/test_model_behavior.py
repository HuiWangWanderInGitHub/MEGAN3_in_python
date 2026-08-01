#!/usr/bin/env python3
"""Behavior tests for the single supported MEGAN calculation workflow."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

from src.megan_model import (
    ModelPaths,
    load_inputs,
    load_namelist,
    precompute_meteorology,
    simulate,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]


class ModelBehaviorTests(unittest.TestCase):
    """Verify maintained indexing, missing-value handling, and switches."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.configuration = load_namelist(PROJECT_ROOT / "MEGAN_namelist.ini")
        cls.bundle = load_inputs(
            cls.configuration.paths,
            cls.configuration.controls,
        )
        cls.result = simulate(cls.bundle, cls.configuration.model)

    def test_missing_lai_does_not_invalidate_the_following_record(self) -> None:
        """One missing LAI record must not propagate into the next record."""

        lai = self.bundle.meteorology["lai"].to_numpy(dtype=float)
        modeled = self.result.all_species.iloc[:, 2].to_numpy(dtype=float)
        missing_indices = np.flatnonzero(~np.isfinite(lai))
        self.assertGreater(missing_indices.size, 0)

        for index in missing_indices:
            if index + 1 < lai.size and np.isfinite(lai[index + 1]):
                self.assertTrue(np.isfinite(modeled[index + 1]))

    def test_isoprene_output_lai_column_contains_lai(self) -> None:
        """The LAI output column must contain LAI rather than wind speed."""

        expected = self.bundle.meteorology["lai"].to_numpy(dtype=float)
        actual = self.result.isoprene["lai_m2_m2"].to_numpy(dtype=float)
        np.testing.assert_allclose(actual, expected, equal_nan=True)

    def test_daily_extrema_ignore_isolated_missing_records(self) -> None:
        """A single NaN must not make a whole day's extreme value missing."""

        precomputed = precompute_meteorology(
            self.bundle.meteorology,
            self.configuration.model,
        )
        self.assertTrue(np.isfinite(precomputed.daily_max_temperature_k[0]))
        self.assertTrue(np.isfinite(precomputed.daily_min_temperature_k[0]))
        self.assertTrue(np.isfinite(precomputed.daily_max_wind_m_s[0]))

    def test_co2_switch_changes_only_isoprene(self) -> None:
        """GAMCO2_YN is independent and applies only to isoprene."""

        controls = replace(self.bundle.controls, co2_response=True)
        result = simulate(
            replace(self.bundle, controls=controls),
            self.configuration.model,
        )

        baseline = self.result.all_species.iloc[:, 2:].to_numpy(dtype=float)
        with_co2 = result.all_species.iloc[:, 2:].to_numpy(dtype=float)
        common = np.isfinite(baseline) & np.isfinite(with_co2)

        isoprene_mask = common[:, 0] & (baseline[:, 0] > 0.0)
        self.assertTrue(
            np.any(
                np.abs(
                    with_co2[isoprene_mask, 0]
                    - baseline[isoprene_mask, 0]
                )
                > 0.0
            )
        )
        np.testing.assert_allclose(
            with_co2[:, 1:][common[:, 1:]],
            baseline[:, 1:][common[:, 1:]],
            rtol=0.0,
            atol=0.0,
        )

    def test_qv_humidity_input_runs(self) -> None:
        """RH_QV=0 must use water-vapor mixing ratio and pressure."""

        with tempfile.TemporaryDirectory() as temporary_directory:
            temporary_path = Path(temporary_directory)
            meteorology = pd.read_csv(
                PROJECT_ROOT / "inputs" / "1.Met_HourlyData_2012_moflux_Kc.csv"
            )
            temperature_c = meteorology["AirTem(degreeC)"]
            saturation_vapor_pressure_pa = (
                0.6112
                * np.exp(17.67 * temperature_c / (temperature_c + 243.5))
                * 1000.0
            )
            vapor_pressure_pa = (
                saturation_vapor_pressure_pa * meteorology["RH(%)"] / 100.0
            )
            molecular_weight_ratio = 18.016 / 28.97
            mixing_ratio = (
                molecular_weight_ratio
                * vapor_pressure_pa
                / (meteorology["AtmPres(Pa)"] - vapor_pressure_pa)
            )
            meteorology = meteorology.drop(columns=["RH(%)"])
            meteorology["QV(kg/kg)"] = mixing_ratio
            meteorology_path = temporary_path / "meteorology_qv.csv"
            meteorology.to_csv(meteorology_path, index=False)

            paths = ModelPaths(
                meteorology=meteorology_path,
                emission_factors=PROJECT_ROOT / "inputs" / "3.EF_LDF.csv",
                pft_fractions=PROJECT_ROOT / "inputs" / "4.PFT_Fraction.csv",
                output_directory=temporary_path / "outputs",
            )
            controls = replace(
                self.configuration.controls,
                humidity_mode="QV",
            )
            result = simulate(
                load_inputs(paths, controls),
                self.configuration.model,
            )
            modeled = result.all_species.iloc[:, 2].to_numpy(dtype=float)
            self.assertGreater(np.count_nonzero(np.isfinite(modeled)), 0)


if __name__ == "__main__":
    unittest.main()
