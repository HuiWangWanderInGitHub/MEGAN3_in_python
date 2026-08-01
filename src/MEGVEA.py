#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""MEGAN emission-activity response functions.

This module contains the environmental activity-factor algorithms used by the
site-scale MEGAN implementation.  The equations are retained from the original
project, while constant coefficient arrays are now created only once at module
import.  This avoids hundreds of small NumPy allocations during a simulation.

The public function names from the original code (for example ``gamma_a`` and
``gamtld``) are preserved for backward compatibility.
"""

from __future__ import annotations

import numpy as np

# The original coefficient tables contain 20 activity classes.  The supplied
# example emission-factor file uses the first 19 classes.
N_ACTIVITY_CLASSES = 20

# Leaf-age response coefficients.
_A_NEW = np.array(
    [
        0.05, 0.05, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.4, 0.4,
        3.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    ],
    dtype=float,
)
_A_GROWING = np.array(
    [
        0.6, 0.6, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 0.6, 0.6,
        3.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    ],
    dtype=float,
)
_A_MATURE = np.ones(N_ACTIVITY_CLASSES, dtype=float)
_A_OLD = np.array(
    [
        0.9, 0.9, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 0.95, 0.95,
        1.2, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    ],
    dtype=float,
)

# Stress-response coefficient tables.
_AIR_QUALITY_COEFFICIENT = np.array(
    [1, 1, 1, 5, 1, 1, 1, 1, 5, 5, 1, 1, 1, 1, 1, 1, 5, 1, 1, 1],
    dtype=float,
)
_HIGH_TEMPERATURE_COEFFICIENT = _AIR_QUALITY_COEFFICIENT.copy()
_LOW_TEMPERATURE_COEFFICIENT = _AIR_QUALITY_COEFFICIENT.copy()
_HIGH_WIND_COEFFICIENT = np.array(
    [1, 1, 5, 5, 5, 5, 5, 5, 5, 5, 1, 1, 1, 1, 1, 1, 5, 1, 1, 1],
    dtype=float,
)

# Temperature-response coefficients.
_CLEO = np.array(
    [
        2.0, 2.0, 1.83, 1.83, 1.83, 1.83, 1.83, 1.83, 2.37, 2.37,
        1.60, 1.83, 2.0, 2.0, 1.83, 1.83, 1.83, 1.83, 1.60, 1.83,
    ],
    dtype=float,
)
_CT1 = np.array(
    [
        95, 95, 80, 80, 80, 80, 80, 80, 130, 130,
        60, 80, 95, 95, 80, 80, 80, 80, 60, 80,
    ],
    dtype=float,
)
_BETA = np.array(
    [
        0.13, 0.13, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.17, 0.17,
        0.08, 0.10, 0.13, 0.13, 0.10, 0.10, 0.10, 0.10, 0.08, 0.10,
    ],
    dtype=float,
)


def _validate_species_index(species_index: int) -> None:
    """Raise a clear error when a species/activity-class index is invalid."""

    if not 0 <= int(species_index) < N_ACTIVITY_CLASSES:
        raise IndexError(
            f"species_index must be between 0 and {N_ACTIVITY_CLASSES - 1}; "
            f"received {species_index}."
        )


def gamma_a(previous_lai: float, current_lai: float, daily_temperature_k: float) -> np.ndarray:
    """Return the leaf-age activity factor for all activity classes.

    Parameters
    ----------
    previous_lai, current_lai
        LAI at the previous and current model records.  The original MEGAN
        formulation assumes an LAI update interval of eight days; the calling
        code is responsible for supplying LAI values with the intended temporal
        interpretation.
    daily_temperature_k
        Mean air temperature for the current day in kelvin.

    Notes
    -----
    The equations and the fixed ``t = 8`` interval are retained from the
    original implementation.  For a zero-LAI canopy, the model driver bypasses
    this function because total emissions are zero.
    """

    lai_previous = float(previous_lai)
    lai_current = float(current_lai)
    temperature = float(daily_temperature_k)
    update_interval_days = 8.0

    if lai_current <= 0.0:
        # The activity factor is immaterial when LAI is zero because emissions
        # are multiplied by LAI.  Returning ones prevents division by zero.
        return np.ones(N_ACTIVITY_CLASSES, dtype=float)

    if lai_previous < lai_current:
        leaf_expansion_days = (
            5.0 + 0.7 * (300.0 - temperature)
            if temperature <= 303.0
            else 2.9
        )
        maturation_days = 2.3 * leaf_expansion_days

        if leaf_expansion_days >= update_interval_days:
            fraction_new = 1.0 - lai_previous / lai_current
        else:
            fraction_new = (
                leaf_expansion_days
                / update_interval_days
                * (1.0 - lai_previous / lai_current)
            )

        if maturation_days >= update_interval_days:
            fraction_mature = lai_previous / lai_current
        else:
            fraction_mature = (
                lai_previous / lai_current
                + (update_interval_days - maturation_days)
                / update_interval_days
                * (1.0 - lai_previous / lai_current)
            )

        fraction_growing = 1.0 - fraction_new - fraction_mature
        fraction_old = 0.0
    elif lai_previous == lai_current:
        fraction_new = 0.0
        fraction_growing = 0.1
        fraction_mature = 0.8
        fraction_old = 0.1
    else:
        fraction_new = 0.0
        fraction_growing = 0.0
        if lai_previous <= 0.0:
            fraction_old = 0.0
            fraction_mature = 1.0
        else:
            fraction_old = (lai_previous - lai_current) / lai_previous
            fraction_mature = 1.0 - fraction_old

    return (
        fraction_new * _A_NEW
        + fraction_growing * _A_GROWING
        + fraction_mature * _A_MATURE
        + fraction_old * _A_OLD
    )


def gamma_cd(layers: int, lai: float) -> np.ndarray:
    """Return the canopy-depth activity factor for each canopy layer."""

    if layers <= 0:
        raise ValueError("layers must be a positive integer")

    if layers == 5:
        canopy_depth = np.array(
            [0.0469101, 0.2307534, 0.5, 0.7692465, 0.9530899],
            dtype=float,
        )
    else:
        canopy_depth = (np.arange(layers, dtype=float) + 0.5) / layers

    lai_depth = np.minimum(float(lai) * canopy_depth, 3.0)
    return -0.2 * lai_depth + 1.3


def gamma_laibidir(lai: float) -> float:
    """Return the bidirectional-exchange LAI activity factor."""

    if lai < 2.0:
        return 0.5 * lai
    if lai <= 6.0:
        return 1.0 - 0.0625 * (lai - 2.0)
    return 0.75


def _piecewise_stress(
    value: float,
    species_index: int,
    threshold: float,
    transition_width: float,
    coefficients: np.ndarray,
    *,
    increasing: bool,
) -> float:
    """Evaluate the common piecewise-linear stress-response shape."""

    _validate_species_index(species_index)
    coefficient = float(coefficients[species_index])

    if increasing:
        if value <= threshold:
            return 1.0
        if value < threshold + transition_width:
            return 1.0 + (coefficient - 1.0) * (
                value - threshold
            ) / transition_width
        return coefficient

    if value >= threshold:
        return 1.0
    if value > threshold - transition_width:
        return 1.0 + (coefficient - 1.0) * (
            threshold - value
        ) / transition_width
    return coefficient


def gamma_aq(species_index: int, air_quality_index: float) -> float:
    """Return the air-quality stress activity factor."""

    return _piecewise_stress(
        air_quality_index,
        species_index,
        threshold=20.0,
        transition_width=30.0,
        coefficients=_AIR_QUALITY_COEFFICIENT,
        increasing=True,
    )


def gamma_ht(species_index: int, maximum_temperature_k: float) -> float:
    """Return the high-temperature stress activity factor."""

    return _piecewise_stress(
        maximum_temperature_k,
        species_index,
        threshold=313.15,
        transition_width=8.0,
        coefficients=_HIGH_TEMPERATURE_COEFFICIENT,
        increasing=True,
    )


def gamma_lt(species_index: int, minimum_temperature_k: float) -> float:
    """Return the low-temperature stress activity factor."""

    return _piecewise_stress(
        minimum_temperature_k,
        species_index,
        threshold=283.15,
        transition_width=8.0,
        coefficients=_LOW_TEMPERATURE_COEFFICIENT,
        increasing=False,
    )


def gamma_hw(species_index: int, maximum_wind_speed_m_s: float) -> float:
    """Return the high-wind stress activity factor."""

    return _piecewise_stress(
        maximum_wind_speed_m_s,
        species_index,
        threshold=12.0,
        transition_width=8.0,
        coefficients=_HIGH_WIND_COEFFICIENT,
        increasing=True,
    )


def gamma_co2(co2_ppm: float) -> float:
    """Return the CO2 inhibition factor for isoprene emissions."""

    isoprene_maximum = 1.344
    co2_exponent = 1.4614
    c_star = 585.0

    if co2_ppm == 400.0:
        return 1.0

    internal_co2 = 0.7 * co2_ppm
    return isoprene_maximum - (
        isoprene_maximum * internal_co2**co2_exponent
    ) / (c_star**co2_exponent + internal_co2**co2_exponent)


def gamtld(
    leaf_temperature_k: float,
    daily_temperature_k: float,
    ten_day_temperature_k: float,
    species_index: int,
) -> float:
    """Light-dependent temperature activity factor."""

    _validate_species_index(species_index)
    if leaf_temperature_k < 260.0:
        return 0.0

    optimum_temperature = 312.5 + 0.6 * (ten_day_temperature_k - 297.0)
    x_value = (1.0 / optimum_temperature - 1.0 / leaf_temperature_k) / 0.00831
    optimum_emission = (
        _CLEO[species_index]
        * np.exp(0.05 * (daily_temperature_k - 297.0))
        * np.exp(0.05 * (ten_day_temperature_k - 297.0))
    )
    ct2 = 230.0
    return float(
        optimum_emission
        * ct2
        * np.exp(_CT1[species_index] * x_value)
        / (
            ct2
            - _CT1[species_index]
            * (1.0 - np.exp(ct2 * x_value))
        )
    )


def gamp(leaf_ppfd: float, daily_mean_ppfd: float) -> float:
    """Light-response activity factor."""

    if daily_mean_ppfd < 0.01:
        return 0.0

    c1 = 1.03
    alpha = 0.004
    return float(
        alpha
        * c1
        * leaf_ppfd
        / np.sqrt(1.0 + alpha**2 * leaf_ppfd**2)
    )


def gamtli(leaf_temperature_k: float, species_index: int) -> float:
    """Light-independent temperature activity factor."""

    _validate_species_index(species_index)
    standard_temperature_k = 303.15
    return float(
        np.exp(
            _BETA[species_index]
            * (leaf_temperature_k - standard_temperature_k)
        )
    )


def ResRC(ppfd: float) -> float:
    """Legacy helper for the HCHO resistance parameterization."""

    return float(
        (0.0027 * 1.066 * ppfd)
        / np.sqrt(1.0 + 0.0027 * 0.0027 * ppfd**2)
    )


def HCHOLeafRes(stomatal_conductance: float) -> float:
    """Legacy HCHO leaf resistance (s cm-1)."""

    return float(0.32 + 50.8 * np.exp(-(stomatal_conductance - 7.5) / 18.884))


def gamhcho(temperature_k: float, ppfd: float, hcho: float) -> tuple[float, int, float]:
    """Return HCHO resistance, exchange direction, and compensation point."""

    stomatal_conductance = ResRC(ppfd) * 65.0
    resistance = HCHOLeafRes(stomatal_conductance)
    delta_temperature = temperature_k - 303.15
    compensation_point = (
        0.552
        + 0.0368 * delta_temperature
        + 0.00363 * delta_temperature**2
        + 0.000189 * delta_temperature**3
    )
    direction = -1 if hcho > compensation_point else 1
    return resistance, direction, float(compensation_point)


def gamma_sm(soil_moisture: float, wilting_point: float) -> float:
    """MEGAN2.1 soil-moisture activity factor."""

    transition_width = 0.04
    upper_threshold = wilting_point + transition_width
    if soil_moisture > upper_threshold:
        return 1.0
    if soil_moisture > wilting_point:
        return float((soil_moisture - wilting_point) / transition_width)
    return 0.0


def gamma_sm_kc(kc: float, kc_max: float, kc_min: float) -> float:
    """ET/PET-based drought activity factor from the project implementation."""

    if not np.isfinite(kc):
        return np.nan
    if kc_max <= kc_min:
        raise ValueError("kc_max must be greater than kc_min")

    k1 = -7.45
    b1 = 3.26
    k2 = -28.76
    b2 = 2.35e6
    maximum_factor = 1.4

    clipped_kc = min(float(kc), float(kc_max))
    normalized_kc = (clipped_kc - kc_min) / (kc_max - kc_min)
    return float(
        maximum_factor
        * (1.0 / (1.0 + b1 * np.exp(k1 * (normalized_kc - 0.2))))
        * (
            (1.0 - 1.0 / maximum_factor)
            / (1.0 + b2 * np.exp(k2 * (1.3 - normalized_kc)))
            + 1.0 / maximum_factor
        )
    )
