#!/usr/bin/env python3
"""Smoke tests covering all six canopy types and day/night branches."""

from __future__ import annotations

import unittest

import numpy as np

from src import MEGCAN as canopy
from src import TIMEFUNC as time_utils


class CanopySmokeTests(unittest.TestCase):
    def test_all_canopy_types_produce_finite_profiles(self) -> None:
        layers = 5
        layer_positions = canopy.gaussian_layer_positions(layers)
        vapor_pressure = canopy.relative_humidity_to_vapor_pressure(60.0, 298.15)

        for canopy_type in range(canopy.N_CANOPY_TYPES):
            for solar, sin_elevation in ((0.0, -0.2), (500.0, 0.7)):
                diffuse_vis, beam_vis, diffuse_nir, beam_nir = (
                    time_utils.partition_solar_radiation(
                        solar,
                        max(1.0, sin_elevation * 1361.5),
                    )
                )
                radiation = canopy.canopy_radiation(
                    layer_positions,
                    layers,
                    3.2,
                    sin_elevation,
                    beam_vis,
                    diffuse_vis,
                    beam_nir,
                    diffuse_nir,
                    canopy_type,
                )
                energy = canopy.canopy_energy_balance(
                    canopy.canopy_temperature_lapse_rate(canopy_type, solar),
                    layers,
                    layer_positions,
                    canopy_type,
                    298.15,
                    2.0,
                    radiation.sun_ppfd,
                    radiation.shade_ppfd,
                    radiation.sun_visible,
                    radiation.shade_visible,
                    radiation.sun_nir,
                    radiation.shade_nir,
                    vapor_pressure,
                )

                for profile in (
                    radiation.sun_fraction,
                    radiation.sun_ppfd,
                    radiation.shade_ppfd,
                    energy.sun_leaf_temperature_k,
                    energy.shade_leaf_temperature_k,
                    energy.wind_speed_m_s,
                ):
                    self.assertEqual(profile.shape, (layers,))
                    self.assertTrue(np.all(np.isfinite(profile)))


if __name__ == "__main__":
    unittest.main()
