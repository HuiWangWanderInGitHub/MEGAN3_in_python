#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Solar-geometry and time-aggregation utilities for site-scale MEGAN.

The original module used Numba decorators on a few small functions.  In this
refactor, the expensive repeated date searches are removed from the model
workflow, so the optional JIT dependency is no longer necessary.  Legacy
function names are retained to avoid breaking existing scripts.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence

import numpy as np

_DEGREES_PER_RADIAN = 57.29578
_PI_APPROX = 3.14159


def solar_elevation_angle(day_of_year: float, latitude_deg: float, hour: float) -> float:
    """Calculate solar elevation angle in degrees.

    Notes
    -----
    The equation is retained from the original MEGAN canopy implementation.
    Despite the old comment calling this a zenith angle, the returned value is
    the *solar elevation* angle because its sine is used directly downstream.
    """

    sin_declination = -np.sin(0.40907) * np.cos(
        6.28 * (day_of_year + 10.0) / 365.0
    )
    cos_declination = np.sqrt(1.0 - sin_declination**2)

    latitude_radians = latitude_deg / _DEGREES_PER_RADIAN
    term_a = np.sin(latitude_radians) * sin_declination
    term_b = np.cos(latitude_radians) * cos_declination
    sin_elevation = term_a + term_b * np.cos(
        2.0 * _PI_APPROX * (hour - 12.0) / 24.0
    )

    # Numerical roundoff can produce values just outside [-1, 1].
    sin_elevation = np.clip(sin_elevation, -1.0, 1.0)
    return float(np.arcsin(sin_elevation) * _DEGREES_PER_RADIAN)


def solar_eccentricity_factor(day_of_year: float) -> float:
    """Return the Earth-Sun distance correction used by the canopy model."""

    return float(1.0 + 0.033 * np.cos(2.0 * 3.14 * (day_of_year - 10.0) / 365.0))


def partition_solar_radiation(solar_w_m2: float, maximum_solar_w_m2: float) -> tuple[float, float, float, float]:
    """Partition solar radiation into direct/diffuse visible and near-IR terms.

    Parameters
    ----------
    solar_w_m2
        Observed incoming solar radiation.
    maximum_solar_w_m2
        Potential maximum solar radiation calculated from solar geometry.

    Returns
    -------
    q_diffuse_visible, q_beam_visible, q_diffuse_nir, q_beam_nir
        Radiation components in W m-2.
    """

    if maximum_solar_w_m2 <= 0.0:
        transmissivity = 0.5
    elif maximum_solar_w_m2 < solar_w_m2:
        transmissivity = 1.0
    else:
        transmissivity = solar_w_m2 / maximum_solar_w_m2

    # Diffuse fraction based on Lizaso et al. (2005), as in the source code.
    diffuse_fraction = 0.156 + 0.86 / (
        1.0 + np.exp(11.1 * (transmissivity - 0.53))
    )

    # Visible (PPFD) fraction based on Goudriaan and van Laar (1994).
    visible_fraction = 0.55 - transmissivity * 0.12

    # Diffuse visible fraction based on Jacovides et al. (2007).
    visible_diffuse_fraction = diffuse_fraction * (
        1.06 + transmissivity * 0.4
    )
    visible_diffuse_fraction = min(visible_diffuse_fraction, 1.0)

    visible = visible_fraction * solar_w_m2
    q_diffuse_visible = visible * visible_diffuse_fraction
    q_beam_visible = visible - q_diffuse_visible

    near_ir = solar_w_m2 - visible
    q_diffuse_nir = near_ir * diffuse_fraction
    q_beam_nir = near_ir - q_diffuse_nir

    return (
        float(q_diffuse_visible),
        float(q_beam_visible),
        float(q_diffuse_nir),
        float(q_beam_nir),
    )


def join_l(values: Iterable[object], separator: str) -> str:
    """Backward-compatible string join helper."""

    return separator.join(str(value) for value in values)


def previous_daily_mean(
    start_day: float,
    current_day: float,
    days: Sequence[float] | np.ndarray,
    values: Sequence[float] | np.ndarray,
    window_days: int,
) -> float:
    """Average prior complete-day means using the original window semantics.

    For the first day, the current day's mean is returned.  For later days, the
    current day is excluded and up to ``window_days`` preceding calendar days
    are averaged.  This reproduces the behavior of the original ``T10`` and
    ``SWC30D`` routines while keeping the logic explicit.
    """

    day_array = np.asarray(days, dtype=float)
    value_array = np.asarray(values, dtype=float)
    day_range = current_day - start_day

    if day_range == 0.0:
        selected = value_array[day_array == current_day]
        return float(np.nanmean(selected))

    first_window_day = (
        current_day - window_days
        if day_range >= float(window_days)
        else start_day
    )
    window = np.arange(first_window_day, current_day)
    daily_means = np.full(window.size, np.nan, dtype=float)

    for index, day in enumerate(window):
        selected = value_array[day_array == day]
        if selected.size:
            daily_means[index] = np.nanmean(selected)

    return float(np.nanmean(daily_means))


def map_previous_daily_mean(
    days: Sequence[float] | np.ndarray,
    values: Sequence[float] | np.ndarray,
    window_days: int,
) -> np.ndarray:
    """Vectorized-by-day wrapper for :func:`previous_daily_mean`.

    Only one rolling value is computed per unique day and then mapped back to
    every record for that day.  This is much faster than repeating the same
    date search for every half-hourly record.
    """

    day_array = np.asarray(days, dtype=float)
    value_array = np.asarray(values, dtype=float)
    if day_array.size != value_array.size:
        raise ValueError("days and values must have the same length")
    if day_array.size == 0:
        return np.array([], dtype=float)

    start_day = float(day_array[0])
    output = np.full(day_array.shape, np.nan, dtype=float)
    for day in np.unique(day_array):
        output[day_array == day] = previous_daily_mean(
            start_day,
            float(day),
            day_array,
            value_array,
            window_days,
        )
    return output


# ---------------------------------------------------------------------------
# Backward-compatible names used by the original main_program.py
# ---------------------------------------------------------------------------

def Calcbeta(Day: float, Lat: float, Hour: float) -> float:
    """Legacy alias of :func:`solar_elevation_angle`."""

    return solar_elevation_angle(Day, Lat, Hour)


def CalcEccentricity(Day: float) -> float:
    """Legacy alias of :func:`solar_eccentricity_factor`."""

    return solar_eccentricity_factor(Day)


def SolarFractions(Solar: float, Maxsolar: float) -> list[float]:
    """Legacy list-returning wrapper for solar-radiation partitioning."""

    return list(partition_solar_radiation(Solar, Maxsolar))


def T10(STIME: float, CTIME: float, TIME: np.ndarray, TEM: np.ndarray) -> float:
    """Legacy 10-day temperature-average wrapper."""

    return previous_daily_mean(STIME, CTIME, TIME, TEM, window_days=10)


def SWC30D(STIME: float, CTIME: float, TIME: np.ndarray, SWC: np.ndarray) -> float:
    """Legacy 30-day soil-water-average wrapper."""

    return previous_daily_mean(STIME, CTIME, TIME, SWC, window_days=30)
