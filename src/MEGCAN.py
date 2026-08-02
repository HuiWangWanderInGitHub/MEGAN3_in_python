#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""MEGAN canopy radiation and leaf-energy-balance model.

The model follows the original site-scale MEGAN canopy formulation, with
three targeted physical corrections: a bounded exponential within-canopy wind
profile, leaf width as the forced-convection length scale, and a temperature
gradient sign consistent with depth measured downward from the canopy top.
The public wrappers are retained so existing scripts remain usable.
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np

# ---------------------------------------------------------------------------
# Canopy characteristics
# ---------------------------------------------------------------------------
# Rows correspond to the 17 characteristics documented below; columns are the
# six canopy types:
#   0 Needleleaf trees
#   1 Tropical forest trees
#   2 Temperate broadleaf trees
#   3 Shrubs
#   4 Herbaceous vegetation
#   5 Crops
#
# Characteristic rows:
#   0 canopy depth (m)
#   1 leaf width (m)
#   2 leaf length (m)
#   3 canopy height (m)
#   4 visible-radiation scattering coefficient
#   5 near-IR scattering coefficient
#   6 diffuse visible reflection coefficient
#   7 diffuse near-IR reflection coefficient
#   8 leaf-clustering coefficient
#   9 leaf IR emissivity
#  10 stomata/cuticle factor
#  11 daytime temperature lapse rate (K m-1)
#  12 nighttime temperature lapse rate (K m-1)
#  13 warm-canopy total humidity change (Pa)
#  14 cool-canopy total humidity change (Pa)
#  15 normalized canopy depth where wind becomes negligible
#  16 canopy transparency
CANOPY_CHARACTERISTICS = np.array(
    [
        [16.0, 16.0, 16.0, 1.0, 0.5, 1.0],
        [0.005, 0.05, 0.05, 0.015, 0.01, 0.02],
        [0.1, 0.1, 0.1, 0.1, 0.15, 0.15],
        [24.0, 24.0, 24.0, 2.0, 0.5, 1.0],
        [0.2, 0.2, 0.2, 0.2, 0.2, 0.2],
        [0.8, 0.8, 0.8, 0.8, 0.8, 0.8],
        [0.057, 0.057, 0.057, 0.057, 0.057, 0.057],
        [0.389, 0.389, 0.389, 0.389, 0.389, 0.389],
        [0.85, 1.1, 0.9, 0.85, 0.7, 0.65],
        [0.95, 0.95, 0.95, 0.95, 0.95, 0.95],
        [1.25, 1.25, 1.25, 1.0, 1.25, 1.25],
        [0.06, 0.06, 0.06, 0.06, 0.06, 0.06],
        [-0.06, -0.06, -0.06, -0.06, -0.06, -0.06],
        [700.0, 700.0, 700.0, 700.0, 700.0, 700.0],
        [150.0, 150.0, 150.0, 150.0, 150.0, 150.0],
        [0.7, 0.7, 0.7, 0.7, 0.7, 0.7],
        [0.2, 0.2, 0.2, 0.2, 0.2, 0.2],
    ],
    dtype=float,
)

N_CANOPY_TYPES = 6
N_CANOPY_CHARACTERISTICS = 17
STEFAN_BOLTZMANN = 5.67e-8  # W m-2 K-4

# Original public constants retained for compatibility.
Canopychar = CANOPY_CHARACTERISTICS
NrTyp = N_CANOPY_TYPES
NrCha = N_CANOPY_CHARACTERISTICS


class CanopyRadiationProfile(NamedTuple):
    """Radiation and sunlit-fraction profiles for one canopy type."""

    sun_fraction: np.ndarray
    absorbed_beam_visible: float
    absorbed_beam_nir: float
    sun_nir: np.ndarray
    shade_nir: np.ndarray
    sun_visible: np.ndarray
    shade_visible: np.ndarray
    sun_ppfd: np.ndarray
    shade_ppfd: np.ndarray
    absorbed_diffuse_visible: np.ndarray
    absorbed_scattered_visible: np.ndarray
    absorbed_diffuse_nir: np.ndarray
    absorbed_scattered_nir: np.ndarray


class LeafEnergyBalanceResult(NamedTuple):
    """Energy-balance solution for a single leaf."""

    outgoing_ir: float
    leaf_temperature_k: float
    sensible_heat: float
    latent_heat: float
    stomatal_resistance: float


class CanopyEnergyBalanceProfile(NamedTuple):
    """Vertical canopy energy-balance profiles."""

    water_vapor_pressure_pa: np.ndarray
    wind_speed_m_s: np.ndarray
    sun_leaf_temperature_k: np.ndarray
    sun_sensible_heat: np.ndarray
    sun_latent_heat: np.ndarray
    sun_net_ir: np.ndarray
    canopy_air_temperature_k: np.ndarray
    shade_leaf_temperature_k: np.ndarray
    shade_sensible_heat: np.ndarray
    shade_latent_heat: np.ndarray
    shade_net_ir: np.ndarray
    sun_stomatal_resistance: np.ndarray
    shade_stomatal_resistance: np.ndarray


def gaussian_layer_positions(number_of_layers: int) -> np.ndarray:
    """Return normalized canopy-depth positions for numerical quadrature."""

    if number_of_layers <= 0:
        raise ValueError("number_of_layers must be positive")
    if number_of_layers == 1:
        return np.array([0.5], dtype=float)
    if number_of_layers == 3:
        return np.array([0.112702, 0.5, 0.887298], dtype=float)
    if number_of_layers == 5:
        return np.array(
            [0.0469101, 0.2307534, 0.5, 0.7692465, 0.9530899],
            dtype=float,
        )
    return (np.arange(number_of_layers, dtype=float) + 0.5) / number_of_layers


def canopy_temperature_lapse_rate(canopy_type: int, solar_w_m2: float) -> float:
    """Return the canopy air-temperature lapse rate for current radiation."""

    transition_solar = 500.0
    daytime_rate = CANOPY_CHARACTERISTICS[11, canopy_type]
    nighttime_rate = CANOPY_CHARACTERISTICS[12, canopy_type]

    if solar_w_m2 > transition_solar:
        return float(daytime_rate)
    if solar_w_m2 > 0.0:
        return float(
            daytime_rate
            - (transition_solar - solar_w_m2)
            / transition_solar
            * (daytime_rate - nighttime_rate)
        )
    return float(nighttime_rate)


def mixing_ratio_to_vapor_pressure(mixing_ratio_kg_kg: float, pressure_pa: float) -> float:
    """Convert water-vapor mixing ratio (kg kg-1) to vapor pressure (Pa)."""

    molecular_weight_ratio = 18.016 / 28.97
    return float(
        mixing_ratio_kg_kg
        / (mixing_ratio_kg_kg + molecular_weight_ratio)
        * pressure_pa
    )


def relative_humidity_to_vapor_pressure(relative_humidity_percent: float, temperature_k: float) -> float:
    """Convert relative humidity (%) to water-vapor pressure (Pa)."""

    temperature_c = temperature_k - 273.15
    saturation_vapor_pressure_kpa = 0.6112 * np.exp(
        17.67 * temperature_c / (temperature_c + 243.5)
    )
    return float(
        saturation_vapor_pressure_kpa
        * relative_humidity_percent
        * 0.01
        * 1000.0
    )


def vapor_pressure_to_relative_humidity(vapor_pressure_kpa: float, temperature_k: float) -> float:
    """Convert water-vapor pressure (kPa) to relative humidity as a fraction."""

    temperature_c = temperature_k - 273.15
    saturation_vapor_pressure_kpa = 0.6112 * np.exp(
        17.67 * temperature_c / (temperature_c + 243.5)
    )
    return float(vapor_pressure_kpa / saturation_vapor_pressure_kpa)


def radiation_extinction_coefficients(
    beam_radiation: float,
    scattering_coefficient: float,
    beam_extinction_black_leaf: float,
    diffuse_extinction_black_leaf: float,
) -> tuple[float, float, float, float]:
    """Calculate beam absorption, reflection, and effective extinction terms."""

    p_value = np.sqrt(1.0 - scattering_coefficient)
    beam_reflection = 1.0 - np.exp(
        (
            -2.0
            * ((1.0 - p_value) / (1.0 + p_value))
            * beam_extinction_black_leaf
        )
        / (1.0 + beam_extinction_black_leaf)
    )

    effective_beam_extinction = beam_extinction_black_leaf * p_value
    effective_diffuse_extinction = diffuse_extinction_black_leaf * p_value
    absorbed_beam = (
        beam_extinction_black_leaf
        * beam_radiation
        * (1.0 - scattering_coefficient)
    )
    return (
        float(absorbed_beam),
        float(beam_reflection),
        float(effective_beam_extinction),
        float(effective_diffuse_extinction),
    )


def absorbed_radiation_components(
    diffuse_radiation: float,
    beam_radiation: float,
    effective_diffuse_extinction: float,
    effective_beam_extinction: float,
    beam_extinction_black_leaf: float,
    scattering_coefficient: float,
    diffuse_reflection: float,
    beam_reflection: float,
    lai_depth: float,
) -> tuple[float, float]:
    """Calculate absorbed diffuse and scattered-beam radiation at LAI depth."""

    absorbed_diffuse = (
        diffuse_radiation
        * effective_diffuse_extinction
        * (1.0 - diffuse_reflection)
        * np.exp(-effective_diffuse_extinction * lai_depth)
    )
    absorbed_scattered = beam_radiation * (
        effective_beam_extinction
        * (1.0 - beam_reflection)
        * np.exp(-effective_beam_extinction * lai_depth)
        - beam_extinction_black_leaf
        * (1.0 - scattering_coefficient)
        * np.exp(-beam_extinction_black_leaf * lai_depth)
    )
    return float(absorbed_diffuse), float(absorbed_scattered)


def canopy_radiation(
    layer_positions: np.ndarray,
    number_of_layers: int,
    lai: float,
    sin_solar_elevation: float,
    beam_visible: float,
    diffuse_visible: float,
    beam_nir: float,
    diffuse_nir: float,
    canopy_type: int,
) -> CanopyRadiationProfile:
    """Calculate visible/near-IR radiation for sunlit and shaded leaves."""

    if number_of_layers != len(layer_positions):
        raise ValueError("number_of_layers must match layer_positions length")
    if not 0 <= canopy_type < N_CANOPY_TYPES:
        raise IndexError(f"canopy_type must be between 0 and {N_CANOPY_TYPES - 1}")

    convert_shade_ppfd = 4.6
    convert_sun_ppfd = 4.0

    sun_fraction = np.zeros(number_of_layers, dtype=float)
    sun_nir = np.zeros(number_of_layers, dtype=float)
    shade_nir = np.zeros(number_of_layers, dtype=float)
    sun_visible = np.zeros(number_of_layers, dtype=float)
    shade_visible = np.zeros(number_of_layers, dtype=float)
    sun_ppfd = np.zeros(number_of_layers, dtype=float)
    shade_ppfd = np.zeros(number_of_layers, dtype=float)
    absorbed_diffuse_visible = np.zeros(number_of_layers, dtype=float)
    absorbed_scattered_visible = np.zeros(number_of_layers, dtype=float)
    absorbed_diffuse_nir = np.zeros(number_of_layers, dtype=float)
    absorbed_scattered_nir = np.zeros(number_of_layers, dtype=float)

    canopy_transparency = CANOPY_CHARACTERISTICS[16, canopy_type]
    adjusted_lai = lai / (1.0 - canopy_transparency)

    is_daytime = (
        beam_visible + diffuse_visible > 0.001
        and sin_solar_elevation > 0.002
        and adjusted_lai > 0.001
    )

    if not is_daytime:
        # Retained from the source code.  The sun/shade distinction is still
        # used by the nighttime IR energy-balance calculation.
        sun_fraction.fill(0.2)
        return CanopyRadiationProfile(
            sun_fraction,
            0.0,
            0.0,
            sun_nir,
            shade_nir,
            sun_visible,
            shade_visible,
            sun_ppfd,
            shade_ppfd,
            absorbed_diffuse_visible,
            absorbed_scattered_visible,
            absorbed_diffuse_nir,
            absorbed_scattered_nir,
        )

    scattering_visible = CANOPY_CHARACTERISTICS[4, canopy_type]
    scattering_nir = CANOPY_CHARACTERISTICS[5, canopy_type]
    diffuse_reflection_visible = CANOPY_CHARACTERISTICS[6, canopy_type]
    diffuse_reflection_nir = CANOPY_CHARACTERISTICS[7, canopy_type]
    clustering = CANOPY_CHARACTERISTICS[8, canopy_type]

    # Black-leaf extinction coefficients, assuming a spherical leaf-angle
    # distribution (0.5 = cos(60 degrees)).
    beam_extinction = clustering * 0.5 / sin_solar_elevation
    diffuse_extinction = 0.8 * clustering

    (
        absorbed_beam_visible,
        beam_reflection_visible,
        effective_beam_visible,
        effective_diffuse_visible,
    ) = radiation_extinction_coefficients(
        beam_visible,
        scattering_visible,
        beam_extinction,
        diffuse_extinction,
    )
    (
        absorbed_beam_nir,
        beam_reflection_nir,
        effective_beam_nir,
        effective_diffuse_nir,
    ) = radiation_extinction_coefficients(
        beam_nir,
        scattering_nir,
        beam_extinction,
        diffuse_extinction,
    )

    for layer in range(number_of_layers):
        lai_depth = adjusted_lai * layer_positions[layer]
        sun_fraction[layer] = np.exp(-beam_extinction * lai_depth)

        diffuse_vis, scattered_vis = absorbed_radiation_components(
            diffuse_visible,
            beam_visible,
            effective_diffuse_visible,
            effective_beam_visible,
            beam_extinction,
            scattering_visible,
            diffuse_reflection_visible,
            beam_reflection_visible,
            lai_depth,
        )
        diffuse_nir_layer, scattered_nir_layer = absorbed_radiation_components(
            diffuse_nir,
            beam_nir,
            effective_diffuse_nir,
            effective_beam_nir,
            beam_extinction,
            scattering_nir,
            diffuse_reflection_nir,
            beam_reflection_nir,
            lai_depth,
        )

        shade_ppfd[layer] = (
            (diffuse_vis + scattered_vis)
            * convert_shade_ppfd
            / (1.0 - scattering_visible)
        )
        sun_ppfd[layer] = shade_ppfd[layer] + (
            absorbed_beam_visible
            * convert_sun_ppfd
            / (1.0 - scattering_visible)
        )

        absorbed_diffuse_visible[layer] = diffuse_vis
        absorbed_scattered_visible[layer] = scattered_vis
        absorbed_diffuse_nir[layer] = diffuse_nir_layer
        absorbed_scattered_nir[layer] = scattered_nir_layer

        shade_visible[layer] = diffuse_vis + scattered_vis
        sun_visible[layer] = shade_visible[layer] + absorbed_beam_visible
        shade_nir[layer] = diffuse_nir_layer + scattered_nir_layer
        sun_nir[layer] = shade_nir[layer] + absorbed_beam_nir

    return CanopyRadiationProfile(
        sun_fraction,
        absorbed_beam_visible,
        absorbed_beam_nir,
        sun_nir,
        shade_nir,
        sun_visible,
        shade_visible,
        sun_ppfd,
        shade_ppfd,
        absorbed_diffuse_visible,
        absorbed_scattered_visible,
        absorbed_diffuse_nir,
        absorbed_scattered_nir,
    )


def stomatal_resistance(ppfd: float) -> float:
    """Leaf stomatal resistance (s m-1)."""

    adjustment = (
        0.0027
        * 1.066
        * ppfd
        / np.sqrt(1.0 + 0.0027 * 0.0027 * ppfd**2)
    )
    return 2000.0 if adjustment < 0.1 else float(200.0 / adjustment)


def latent_heat_of_vaporization(temperature_k: float) -> float:
    """Latent heat of vaporization (J kg-1)."""

    return float(2_501_000.0 - 2_370.0 * (temperature_k - 273.15))


def leaf_ir_emission(temperature_k: float, emissivity: float) -> float:
    """Two-sided leaf infrared emission (W m-2)."""

    return float(emissivity * STEFAN_BOLTZMANN * 2.0 * temperature_k**4)


def exposed_leaf_ir_input(vapor_pressure_pa: float, temperature_k: float) -> float:
    """Downward atmospheric IR incident on an exposed leaf surface."""

    atmospheric_emissivity = (
        0.7
        + 5.95
        * (vapor_pressure_pa / 1000.0)
        * 1.0e-4
        * np.exp(1500.0 / temperature_k)
    )
    return float(atmospheric_emissivity * STEFAN_BOLTZMANN * temperature_k**4)


def leaf_boundary_layer_conductance(
    forced_conductance: float,
    leaf_air_temperature_difference: float,
    leaf_length_m: float,
) -> float:
    """Combined forced and free-convection heat-transfer coefficient."""

    if leaf_air_temperature_difference >= 0.0:
        free_conductance = (
            0.5
            * 0.00253
            * (
                160_000_000.0
                * leaf_air_temperature_difference
                / leaf_length_m**3
            )
            ** 0.25
            / leaf_length_m
        )
    else:
        free_conductance = 0.0
    return float(forced_conductance + free_conductance)


def leaf_sensible_heat(
    leaf_air_temperature_difference: float,
    conductance: float,
) -> float:
    """Two-sided convective sensible heat flux (W m-2)."""

    return float(2.0 * conductance * leaf_air_temperature_difference)


def leaf_latent_heat(
    leaf_temperature_k: float,
    ambient_vapor_density_kg_m3: float,
    latent_heat_j_kg: float,
    heat_conductance: float,
    stomatal_resistance_s_m: float,
    transpiration_type: float,
) -> float:
    """Latent heat term in the leaf energy balance (W m-2)."""

    leaf_resistance = (
        1.0 / (1.075 * (heat_conductance / 1231.0))
        + stomatal_resistance_s_m
    )
    saturation_vapor_pressure = 10.0 ** (
        -2937.4 / leaf_temperature_k
        - 4.9283 * np.log10(leaf_temperature_k)
        + 23.5518
    )
    saturation_vapor_density = (
        0.2165 * saturation_vapor_pressure / leaf_temperature_k
    )
    vapor_deficit = saturation_vapor_density - ambient_vapor_density_kg_m3
    latent_heat = (
        transpiration_type
        / leaf_resistance
        * latent_heat_j_kg
        * vapor_deficit
    )
    return float(max(latent_heat, 0.0))


def leaf_energy_balance(
    ppfd: float,
    absorbed_shortwave: float,
    incoming_ir: float,
    emissivity: float,
    transpiration_type: float,
    leaf_width_m: float,
    leaf_length_m: float,
    air_temperature_k: float,
    vapor_pressure_pa: float,
    wind_speed_m_s: float,
) -> LeafEnergyBalanceResult:
    """Solve leaf temperature and energy fluxes by iterative energy balance."""

    effective_wind = max(float(wind_speed_m_s), 0.001)
    ambient_vapor_density = 0.002165 * vapor_pressure_pa / air_temperature_k

    if leaf_width_m <= 0.0 or leaf_length_m <= 0.0:
        raise ValueError("leaf_width_m and leaf_length_m must be positive")

    # Forced convection is controlled by the leaf dimension normal to the
    # airflow.  Use leaf width here; leaf length remains the characteristic
    # scale for the free-convection term below.  The previous implementation
    # accepted leaf_width_m but never used it.
    forced_conductance = 0.0259 / (
        0.004 * np.sqrt(leaf_width_m / effective_wind)
    )
    resistance = stomatal_resistance(ppfd)
    outgoing_ir_at_air_temperature = leaf_ir_emission(
        air_temperature_k,
        emissivity,
    )
    latent_heat_j_kg = latent_heat_of_vaporization(air_temperature_k)
    latent_heat_at_air_temperature = leaf_latent_heat(
        air_temperature_k,
        ambient_vapor_density,
        latent_heat_j_kg,
        forced_conductance,
        resistance,
        transpiration_type,
    )

    residual_at_air_temperature = (
        absorbed_shortwave
        + incoming_ir
        - outgoing_ir_at_air_temperature
        - latent_heat_at_air_temperature
    )
    if residual_at_air_temperature == 0.0:
        residual_at_air_temperature = -1.0

    temperature_difference = 1.0
    balance = 10.0
    for _ in range(10):
        if abs(balance) <= 2.0:
            break

        conductance = leaf_boundary_layer_conductance(
            forced_conductance,
            temperature_difference,
            leaf_length_m,
        )
        sensible_heat = leaf_sensible_heat(temperature_difference, conductance)
        latent_heat_j_kg = latent_heat_of_vaporization(
            air_temperature_k + temperature_difference
        )
        latent_heat = leaf_latent_heat(
            air_temperature_k + temperature_difference,
            ambient_vapor_density,
            latent_heat_j_kg,
            conductance,
            resistance,
            transpiration_type,
        )
        latent_heat_change = latent_heat - latent_heat_at_air_temperature

        outgoing_ir = leaf_ir_emission(
            air_temperature_k + temperature_difference,
            emissivity,
        )
        outgoing_ir_change = outgoing_ir - outgoing_ir_at_air_temperature

        denominator = (
            sensible_heat + latent_heat_change + outgoing_ir_change
        ) / temperature_difference
        temperature_difference = residual_at_air_temperature / denominator
        balance = (
            absorbed_shortwave
            + incoming_ir
            - outgoing_ir
            - sensible_heat
            - latent_heat
        )

    temperature_difference = float(np.clip(temperature_difference, -10.0, 10.0))
    leaf_temperature = air_temperature_k + temperature_difference
    conductance = leaf_boundary_layer_conductance(
        forced_conductance,
        temperature_difference,
        leaf_length_m,
    )
    sensible_heat = leaf_sensible_heat(temperature_difference, conductance)
    latent_heat = leaf_latent_heat(
        leaf_temperature,
        ambient_vapor_density,
        latent_heat_j_kg,
        conductance,
        resistance,
        transpiration_type,
    )
    outgoing_ir = leaf_ir_emission(leaf_temperature, emissivity)

    return LeafEnergyBalanceResult(
        outgoing_ir,
        leaf_temperature,
        sensible_heat,
        latent_heat,
        resistance,
    )


def canopy_energy_balance(
    temperature_lapse_rate: float,
    number_of_layers: int,
    layer_positions: np.ndarray,
    canopy_type: int,
    above_canopy_temperature_k: float,
    above_canopy_wind_m_s: float,
    sun_ppfd: np.ndarray,
    shade_ppfd: np.ndarray,
    sun_visible: np.ndarray,
    shade_visible: np.ndarray,
    sun_nir: np.ndarray,
    shade_nir: np.ndarray,
    above_canopy_vapor_pressure_pa: float,
) -> CanopyEnergyBalanceProfile:
    """Calculate vertical air, wind, humidity, and leaf-temperature profiles."""

    water_vapor_pressure = np.zeros(number_of_layers, dtype=float)
    wind_speed = np.zeros(number_of_layers, dtype=float)
    sun_leaf_temperature = np.zeros(number_of_layers, dtype=float)
    sun_sensible_heat = np.zeros(number_of_layers, dtype=float)
    sun_latent_heat = np.zeros(number_of_layers, dtype=float)
    sun_net_ir = np.zeros(number_of_layers, dtype=float)
    canopy_air_temperature = np.zeros(number_of_layers, dtype=float)
    shade_leaf_temperature = np.zeros(number_of_layers, dtype=float)
    shade_sensible_heat = np.zeros(number_of_layers, dtype=float)
    shade_latent_heat = np.zeros(number_of_layers, dtype=float)
    shade_net_ir = np.zeros(number_of_layers, dtype=float)
    sun_stomatal_resistance = np.zeros(number_of_layers, dtype=float)
    shade_stomatal_resistance = np.zeros(number_of_layers, dtype=float)

    canopy_depth = CANOPY_CHARACTERISTICS[0, canopy_type]
    leaf_width = CANOPY_CHARACTERISTICS[1, canopy_type]
    leaf_length = CANOPY_CHARACTERISTICS[2, canopy_type]
    canopy_height = CANOPY_CHARACTERISTICS[3, canopy_type]
    emissivity = CANOPY_CHARACTERISTICS[9, canopy_type]
    transpiration_type = CANOPY_CHARACTERISTICS[10, canopy_type]
    warm_humidity_change = CANOPY_CHARACTERISTICS[13, canopy_type]
    cool_humidity_change = CANOPY_CHARACTERISTICS[14, canopy_type]
    normalized_no_wind_depth = CANOPY_CHARACTERISTICS[15, canopy_type]

    if above_canopy_temperature_k > 288.0:
        humidity_gradient = warm_humidity_change / canopy_height
    elif above_canopy_temperature_k > 278.0:
        humidity_gradient = (
            warm_humidity_change
            - (288.0 - above_canopy_temperature_k)
            / 10.0
            * (warm_humidity_change - cool_humidity_change)
        ) / canopy_height
    else:
        humidity_gradient = cool_humidity_change / canopy_height

    layer_depth_m = canopy_depth * layer_positions

    # layer_depth_m increases downward from the canopy top.  Positive daytime
    # lapse rates therefore represent cooling with depth, while the negative
    # nighttime rate represents a warmer canopy interior.
    canopy_air_temperature[:] = (
        above_canopy_temperature_k
        - temperature_lapse_rate * layer_depth_m
    )
    water_vapor_pressure[:] = (
        above_canopy_vapor_pressure_pa
        + humidity_gradient * layer_depth_m
    )

    # Use a bounded exponential attenuation with normalized canopy depth.
    # The characteristic in row 15 denotes the depth fraction at which only
    # five percent of the above-floor wind excess remains.  This guarantees a
    # positive, monotonic profile that never exceeds the above-canopy wind.
    top_wind = max(float(above_canopy_wind_m_s), 0.001)
    minimum_wind = min(0.05, top_wind)
    no_wind_depth = max(float(normalized_no_wind_depth), 1.0e-6)
    attenuation_coefficient = -np.log(0.05) / no_wind_depth
    normalized_depth = layer_depth_m / max(canopy_depth, 1.0e-6)
    wind_speed[:] = minimum_wind + (top_wind - minimum_wind) * np.exp(
        -attenuation_coefficient * normalized_depth
    )

    for layer in range(number_of_layers):
        atmospheric_emissivity = (
            0.7
            + 5.95
            * (water_vapor_pressure[layer] / 1000.0)
            * 1.0e-4
            * np.exp(1500.0 / canopy_air_temperature[layer])
        )
        unexposed_ir = leaf_ir_emission(
            canopy_air_temperature[layer],
            atmospheric_emissivity,
        )
        shade_incoming_ir = unexposed_ir
        sun_incoming_ir = (
            0.75 * unexposed_ir
            + 0.5
            * exposed_leaf_ir_input(
                above_canopy_vapor_pressure_pa,
                above_canopy_temperature_k,
            )
        )

        sun_result = leaf_energy_balance(
            sun_ppfd[layer],
            sun_visible[layer] + sun_nir[layer],
            sun_incoming_ir,
            emissivity,
            transpiration_type,
            leaf_width,
            leaf_length,
            canopy_air_temperature[layer],
            water_vapor_pressure[layer],
            wind_speed[layer],
        )
        sun_leaf_temperature[layer] = sun_result.leaf_temperature_k
        sun_sensible_heat[layer] = sun_result.sensible_heat
        sun_latent_heat[layer] = sun_result.latent_heat
        sun_net_ir[layer] = sun_incoming_ir - sun_result.outgoing_ir
        sun_stomatal_resistance[layer] = sun_result.stomatal_resistance

        shade_result = leaf_energy_balance(
            shade_ppfd[layer],
            shade_visible[layer] + shade_nir[layer],
            shade_incoming_ir,
            emissivity,
            transpiration_type,
            leaf_width,
            leaf_length,
            canopy_air_temperature[layer],
            water_vapor_pressure[layer],
            wind_speed[layer],
        )
        shade_leaf_temperature[layer] = shade_result.leaf_temperature_k
        shade_sensible_heat[layer] = shade_result.sensible_heat
        shade_latent_heat[layer] = shade_result.latent_heat
        shade_net_ir[layer] = shade_incoming_ir - shade_result.outgoing_ir
        shade_stomatal_resistance[layer] = shade_result.stomatal_resistance

    return CanopyEnergyBalanceProfile(
        water_vapor_pressure,
        wind_speed,
        sun_leaf_temperature,
        sun_sensible_heat,
        sun_latent_heat,
        sun_net_ir,
        canopy_air_temperature,
        shade_leaf_temperature,
        shade_sensible_heat,
        shade_latent_heat,
        shade_net_ir,
        sun_stomatal_resistance,
        shade_stomatal_resistance,
    )


def adjusted_lai(lai: float, canopy_type: int) -> float:
    """Return LAI adjusted for canopy transparency."""

    return float(lai / (1.0 - CANOPY_CHARACTERISTICS[16, canopy_type]))


# ---------------------------------------------------------------------------
# Backward-compatible wrappers using the original public names and return order
# ---------------------------------------------------------------------------

def GaussianDist(NLayers: int) -> np.ndarray:
    """Legacy alias of :func:`gaussian_layer_positions`."""

    return gaussian_layer_positions(NLayers)


def Stability(Cantype: int, Solar: float) -> float:
    """Legacy alias of :func:`canopy_temperature_lapse_rate`."""

    return canopy_temperature_lapse_rate(Cantype, Solar)


def WaterVapPres(Dens: float, Pres: float) -> float:
    """Legacy alias converting mixing ratio to vapor pressure."""

    return mixing_ratio_to_vapor_pressure(Dens, Pres)


def RHtoWaterVapPres(RH: float, T: float) -> float:
    """Legacy alias converting relative humidity to vapor pressure."""

    return relative_humidity_to_vapor_pressure(RH, T)


def WaterVapPrestoRH(WVP: float, T: float) -> float:
    """Legacy alias converting vapor pressure to relative humidity."""

    return vapor_pressure_to_relative_humidity(WVP, T)


def CalcExtCoeff(Qbeam: float, scat: float, kb: float, kd: float) -> list[float]:
    """Legacy list-returning wrapper for extinction coefficients."""

    return list(radiation_extinction_coefficients(Qbeam, scat, kb, kd))


def CalcRadComponents(
    Qdiff: float,
    Qbeam: float,
    Kdp: float,
    Kbp: float,
    Kb: float,
    Scat: float,
    Refld: float,
    Reflb: float,
    LAIdepth: float,
) -> tuple[float, float]:
    """Legacy wrapper for absorbed diffuse and scattered radiation."""

    return absorbed_radiation_components(
        Qdiff, Qbeam, Kdp, Kbp, Kb, Scat, Refld, Reflb, LAIdepth
    )


def CanopyRad(
    Distgauss: np.ndarray,
    layers: int,
    LAI: float,
    Sinbeta: float,
    Qbeamv: float,
    Qdiffv: float,
    Qbeamn: float,
    Qdiffn: float,
    Cantype: int,
) -> list[object]:
    """Legacy list-returning wrapper for canopy radiation profiles."""

    profile = list(
        canopy_radiation(
            Distgauss,
            layers,
            LAI,
            Sinbeta,
            Qbeamv,
            Qdiffv,
            Qbeamn,
            Qdiffn,
            Cantype,
        )
    )
    # Preserve an original return-type inconsistency for older external code:
    # at night QbAbsN remained its initialized zero array, whereas QbAbsV was
    # replaced by scalar zero.  The modern ``canopy_radiation`` API returns
    # scalars consistently; only this compatibility wrapper reproduces the old
    # nighttime list shape.
    adjusted = LAI / (1.0 - CANOPY_CHARACTERISTICS[16, Cantype])
    is_daytime = (
        Qbeamv + Qdiffv > 0.001
        and Sinbeta > 0.002
        and adjusted > 0.001
    )
    if not is_daytime:
        profile[2] = np.zeros(layers, dtype=float)
    return profile


def ResSC(PPFD: float) -> float:
    """Legacy alias of :func:`stomatal_resistance`."""

    return stomatal_resistance(PPFD)


def LHV(Tk: float) -> float:
    """Legacy alias of :func:`latent_heat_of_vaporization`."""

    return latent_heat_of_vaporization(Tk)


def LeafIR(Tk: float, Eps: float) -> float:
    """Legacy alias of :func:`leaf_ir_emission`."""

    return leaf_ir_emission(Tk, Eps)


def ExposedLeafIRin(HumidPa: float, Tk: float) -> float:
    """Legacy alias of :func:`exposed_leaf_ir_input`."""

    return exposed_leaf_ir_input(HumidPa, Tk)


def LeafBLC(GHforced: float, Tdelta: float, Llength: float) -> float:
    """Legacy alias of :func:`leaf_boundary_layer_conductance`."""

    return leaf_boundary_layer_conductance(GHforced, Tdelta, Llength)


def LeafH(Tdelta: float, GH: float) -> float:
    """Legacy alias of :func:`leaf_sensible_heat`."""

    return leaf_sensible_heat(Tdelta, GH)


def LeafLE(
    Tleaf: float,
    Ambvap: float,
    LatHv: float,
    GH: float,
    StomRes: float,
    TranspireType: float,
) -> float:
    """Legacy alias of :func:`leaf_latent_heat`."""

    return leaf_latent_heat(Tleaf, Ambvap, LatHv, GH, StomRes, TranspireType)


def LeafEB(
    PPFD: float,
    Q: float,
    IRin: float,
    Eps: float,
    TranspireType: float,
    Lwidth: float,
    Llength: float,
    TairK: float,
    HumidairPa: float,
    Ws: float,
) -> list[float]:
    """Legacy list-returning wrapper for one-leaf energy balance."""

    return list(
        leaf_energy_balance(
            PPFD,
            Q,
            IRin,
            Eps,
            TranspireType,
            Lwidth,
            Llength,
            TairK,
            HumidairPa,
            Ws,
        )
    )


def CanopyEB(
    Trate: float,
    Layers: int,
    Distgauss: np.ndarray,
    Cantype: int,
    TairK0: float,
    Ws0: float,
    SunPPFD: np.ndarray,
    ShadePPFD: np.ndarray,
    SunQv: np.ndarray,
    ShadeQv: np.ndarray,
    SunQn: np.ndarray,
    ShadeQn: np.ndarray,
    HumidairPa0: float,
) -> list[np.ndarray]:
    """Legacy list-returning wrapper for canopy energy balance."""

    return list(
        canopy_energy_balance(
            Trate,
            Layers,
            Distgauss,
            Cantype,
            TairK0,
            Ws0,
            SunPPFD,
            ShadePPFD,
            SunQv,
            ShadeQv,
            SunQn,
            ShadeQn,
            HumidairPa0,
        )
    )


def laiadj(laic: float, S: int) -> float:
    """Legacy alias of :func:`adjusted_lai`."""

    return adjusted_lai(laic, S)


def CalcwaterVPpa(RH: float, tk: float) -> float:
    """Legacy Tetens helper retained with its original kPa return units."""

    return float(
        RH
        * 0.01
        * (
            0.6112
            * np.exp((17.67 * (tk - 273.16)) / (tk - 29.66))
        )
    )
