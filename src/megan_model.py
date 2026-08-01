#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""High-level site-scale MEGAN model driver.

This module separates file I/O, validation, meteorological preprocessing,
canopy calculations, emission calculations, diagnostics, and output writing.
The scientific response functions remain in :mod:`src.MEGVEA`, while canopy
physics remain in :mod:`src.MEGCAN`.
"""

from __future__ import annotations

import configparser
from dataclasses import dataclass
import json
import os
from pathlib import Path
import shutil
import warnings
from typing import Literal, NamedTuple

import numpy as np
import pandas as pd

from . import MEGCAN as canopy
from . import MEGVEA as activity
from . import TIMEFUNC as time_utils


# Five-point Gaussian quadrature weights retained from the original model.
_LAYER_WEIGHTS_5 = np.array(
    [0.1184635, 0.2393144, 0.284444444, 0.2393144, 0.1184635],
    dtype=float,
)


@dataclass(frozen=True)
class ModelPaths:
    """Input and output paths for one model run."""

    meteorology: Path
    emission_factors: Path
    pft_fractions: Path
    output_directory: Path


@dataclass(frozen=True)
class ModelOptions:
    """Numerical and diagnostic options that are not stored in CSV files."""

    number_of_layers: int = 5
    solar_constant_w_m2: float = 1361.5
    solar_to_ppfd: float = 2.1
    air_quality_index: float = 40.0
    co2_ppm: float = 403.0
    kc_min: float = 0.0
    kc_max: float = 0.82
    isoprene_molecular_weight_g_mol: float = 68.12
    daytime_start_hour: float = 9.0
    daytime_end_hour: float = 17.0

    def __post_init__(self) -> None:
        """Validate option ranges immediately after dataclass construction."""

        if self.number_of_layers <= 0:
            raise ValueError("number_of_layers must be positive")
        if self.solar_constant_w_m2 <= 0.0:
            raise ValueError("solar_constant_w_m2 must be positive")
        if self.solar_to_ppfd <= 0.0:
            raise ValueError("solar_to_ppfd must be positive")
        if self.air_quality_index < 0.0:
            raise ValueError("air_quality_index cannot be negative")
        if self.kc_max <= self.kc_min:
            raise ValueError("kc_max must be greater than kc_min")
        if self.isoprene_molecular_weight_g_mol <= 0.0:
            raise ValueError("isoprene_molecular_weight_g_mol must be positive")
        if not 0.0 <= self.daytime_start_hour <= 24.0:
            raise ValueError("daytime_start_hour must be between 0 and 24")
        if not 0.0 <= self.daytime_end_hour <= 24.0:
            raise ValueError("daytime_end_hour must be between 0 and 24")
        if self.daytime_end_hour < self.daytime_start_hour:
            raise ValueError("daytime_end_hour must be >= daytime_start_hour")
        if self.co2_ppm <= 0.0:
            raise ValueError("co2_ppm must be positive")


@dataclass(frozen=True)
class OutputOptions:
    """Output filenames, formatting, and diagnostic-figure settings."""

    all_species_filename: str = "Moflux_2012_simulation_all.csv"
    isoprene_filename: str = "Moflux_2012_simulation_isoprene.csv"
    metrics_filename: str = "Moflux_2012_metrics.json"
    diagnostic_figure_filename: str = "Moflux_2012_isoprene_diagnostics.png"
    float_format: str = "%.6f"
    missing_value: str = "nan"
    create_diagnostic_figure: bool = True
    show_diagnostic_figure: bool = False
    copy_namelist_to_output: bool = True
    namelist_copy_filename: str = "MEGAN_namelist_used.ini"

    def __post_init__(self) -> None:
        """Validate output names and formatting settings."""

        filenames = {
            "all_species_filename": self.all_species_filename,
            "isoprene_filename": self.isoprene_filename,
            "metrics_filename": self.metrics_filename,
            "diagnostic_figure_filename": self.diagnostic_figure_filename,
            "namelist_copy_filename": self.namelist_copy_filename,
        }
        for label, filename in filenames.items():
            if not filename.strip():
                raise ValueError(f"{label} cannot be empty")

        active_names = [
            self.all_species_filename,
            self.isoprene_filename,
            self.metrics_filename,
            self.diagnostic_figure_filename,
        ]
        if self.copy_namelist_to_output:
            active_names.append(self.namelist_copy_filename)
        if len(set(active_names)) != len(active_names):
            raise ValueError("Output filenames must be distinct")
        if self.show_diagnostic_figure and not self.create_diagnostic_figure:
            raise ValueError(
                "show_diagnostic_figure requires create_diagnostic_figure=true"
            )
        try:
            self.float_format % 1.23456789
        except (TypeError, ValueError) as error:
            raise ValueError(
                "float_format must be a valid printf-style format, e.g. %.6f"
            ) from error


@dataclass(frozen=True)
class ControlParameters:
    """Activity switches and site parameters for one model run."""

    bidirectional_lai_response: bool
    air_quality_response: bool
    high_temperature_response: bool
    low_temperature_response: bool
    high_wind_response: bool
    soil_moisture_response: bool
    co2_response: bool
    humidity_mode: Literal["RH", "QV"]
    latitude_deg: float
    wilting_point: float


@dataclass(frozen=True)
class RunConfiguration:
    """Complete run configuration loaded from one namelist-style INI file."""

    source_path: Path
    paths: ModelPaths
    controls: ControlParameters
    model: ModelOptions
    output: OutputOptions


@dataclass(frozen=True)
class SpeciesParameters:
    """Species/activity-class metadata read from ``3.EF_LDF.csv``."""

    names: tuple[str, ...]
    emission_factors_nmol_m2_s: np.ndarray
    light_dependent_fraction: np.ndarray


@dataclass(frozen=True)
class PFTParameters:
    """Canopy-type fractions read from ``4.PFT_Fraction.csv``."""

    names: tuple[str, ...]
    fractions_percent: np.ndarray


@dataclass(frozen=True)
class InputBundle:
    """Validated inputs required by the model."""

    meteorology: pd.DataFrame
    controls: ControlParameters
    species: SpeciesParameters
    pfts: PFTParameters


@dataclass(frozen=True)
class PrecomputedMeteorology:
    """Daily and rolling quantities mapped to each input record."""

    daily_mean_ppfd: np.ndarray
    daily_mean_temperature_k: np.ndarray
    previous_10_day_temperature_k: np.ndarray
    daily_max_temperature_k: np.ndarray
    daily_min_temperature_k: np.ndarray
    daily_max_wind_m_s: np.ndarray


class CanopyState(NamedTuple):
    """PFT-weighted sun/shade canopy state for all vertical layers."""

    sun_leaf_temperature_k: np.ndarray
    shade_leaf_temperature_k: np.ndarray
    sun_ppfd: np.ndarray
    shade_ppfd: np.ndarray
    sun_fraction: np.ndarray


@dataclass(frozen=True)
class SimulationResult:
    """Data tables and diagnostics produced by one simulation."""

    all_species: pd.DataFrame
    isoprene: pd.DataFrame
    metrics: dict[str, float | int | str | None]


def _clean_text(value: object) -> str:
    """Strip whitespace and byte-order marks from CSV labels."""

    return str(value).replace("\ufeff", "").strip()


def _namelist_section(
    parser: configparser.ConfigParser,
    section_name: str,
) -> configparser.SectionProxy:
    """Return a namelist section using a case-insensitive section match."""

    for existing_name in parser.sections():
        if existing_name.casefold() == section_name.casefold():
            return parser[existing_name]
    raise ValueError(
        f"Required namelist section '[{section_name}]' is missing"
    )


def _namelist_path(raw_value: str, base_directory: Path) -> Path:
    """Expand a path from the namelist and resolve it relative to that file."""

    expanded = os.path.expandvars(os.path.expanduser(raw_value.strip()))
    if not expanded:
        raise ValueError("A namelist path value cannot be empty")
    path = Path(expanded)
    if not path.is_absolute():
        path = base_directory / path
    return path.resolve()


def _namelist_humidity_mode(raw_value: str) -> Literal["RH", "QV"]:
    """Translate the familiar RH_QV namelist flag to an explicit mode."""

    normalized = raw_value.strip().casefold()
    if normalized in {"1", "rh", "relative_humidity", "relative humidity"}:
        return "RH"
    if normalized in {
        "0",
        "qv",
        "mixing_ratio",
        "water_vapor_mixing_ratio",
        "water vapor mixing ratio",
    }:
        return "QV"
    raise ValueError(
        "SITE.RH_QV must be 1/RH for relative humidity or 0/QV for "
        f"water-vapor mixing ratio; received {raw_value!r}"
    )


def load_namelist(path: Path) -> RunConfiguration:
    """Load all run paths and settings from one namelist-style INI file.

    Relative paths are interpreted relative to the namelist location, not the
    shell's current working directory.  This makes a project directory portable
    and allows ``python main_program.py`` to be called from anywhere.
    """

    source_path = path.expanduser().resolve()
    if not source_path.exists():
        raise FileNotFoundError(f"Namelist file not found: {source_path}")

    parser = configparser.ConfigParser(interpolation=None)
    try:
        with source_path.open("r", encoding="utf-8-sig") as file_handle:
            parser.read_file(file_handle)
    except configparser.Error as error:
        raise ValueError(f"Invalid namelist syntax in {source_path}: {error}") from error

    files = _namelist_section(parser, "FILES")
    model = _namelist_section(parser, "MODEL")
    site = _namelist_section(parser, "SITE")
    switches = _namelist_section(parser, "ACTIVITY_SWITCHES")
    diagnostics = _namelist_section(parser, "DIAGNOSTICS")
    output = _namelist_section(parser, "OUTPUT")
    base_directory = source_path.parent

    def required_text(
        section: configparser.SectionProxy,
        key: str,
    ) -> str:
        """Read a required non-empty string with a targeted error message."""

        if key not in section:
            raise ValueError(
                f"Required namelist item '[{section.name}] {key}' is missing"
            )
        value = section.get(key, raw=True).strip()
        if not value:
            raise ValueError(
                f"Namelist item '[{section.name}] {key}' cannot be empty"
            )
        return value

    def required_float(
        section: configparser.SectionProxy,
        key: str,
    ) -> float:
        """Read a required floating-point value."""

        required_text(section, key)
        try:
            return section.getfloat(key)
        except ValueError as error:
            raise ValueError(
                f"Namelist item '[{section.name}] {key}' must be numeric"
            ) from error

    def required_int(
        section: configparser.SectionProxy,
        key: str,
    ) -> int:
        """Read a required integer value."""

        required_text(section, key)
        try:
            return section.getint(key)
        except ValueError as error:
            raise ValueError(
                f"Namelist item '[{section.name}] {key}' must be an integer"
            ) from error

    def required_bool(
        section: configparser.SectionProxy,
        key: str,
    ) -> bool:
        """Read a required boolean, accepting 0/1 and true/false."""

        required_text(section, key)
        try:
            return section.getboolean(key)
        except ValueError as error:
            raise ValueError(
                f"Namelist item '[{section.name}] {key}' must be 0/1, "
                "true/false, yes/no, or on/off"
            ) from error

    paths = ModelPaths(
        meteorology=_namelist_path(
            required_text(files, "meteorology"), base_directory
        ),
        emission_factors=_namelist_path(
            required_text(files, "emission_factors"), base_directory
        ),
        pft_fractions=_namelist_path(
            required_text(files, "pft_fractions"), base_directory
        ),
        output_directory=_namelist_path(
            required_text(files, "output_directory"), base_directory
        ),
    )

    controls = ControlParameters(
        bidirectional_lai_response=required_bool(switches, "GAMBD_YN"),
        air_quality_response=required_bool(switches, "GAMAQ_YN"),
        high_temperature_response=required_bool(switches, "GAMHT_YN"),
        low_temperature_response=required_bool(switches, "GAMLT_YN"),
        high_wind_response=required_bool(switches, "GAMHW_YN"),
        soil_moisture_response=required_bool(switches, "GAMSM_YN"),
        co2_response=required_bool(switches, "GAMCO2_YN"),
        humidity_mode=_namelist_humidity_mode(required_text(site, "RH_QV")),
        latitude_deg=required_float(site, "Latitude"),
        wilting_point=required_float(site, "WT"),
    )
    if not -90.0 <= controls.latitude_deg <= 90.0:
        raise ValueError(
            "SITE.Latitude must be between -90 and 90 degrees; received "
            f"{controls.latitude_deg}"
        )

    model_options = ModelOptions(
        number_of_layers=required_int(model, "number_of_layers"),
        solar_constant_w_m2=required_float(model, "solar_constant_w_m2"),
        solar_to_ppfd=required_float(model, "solar_to_ppfd"),
        air_quality_index=required_float(model, "air_quality_index"),
        co2_ppm=required_float(model, "co2_ppm"),
        kc_min=required_float(model, "kc_min"),
        kc_max=required_float(model, "kc_max"),
        isoprene_molecular_weight_g_mol=required_float(
            model, "isoprene_molecular_weight_g_mol"
        ),
        daytime_start_hour=required_float(
            diagnostics, "daytime_start_hour"
        ),
        daytime_end_hour=required_float(diagnostics, "daytime_end_hour"),
    )

    output_options = OutputOptions(
        all_species_filename=required_text(output, "all_species_filename"),
        isoprene_filename=required_text(output, "isoprene_filename"),
        metrics_filename=required_text(output, "metrics_filename"),
        diagnostic_figure_filename=required_text(
            output, "diagnostic_figure_filename"
        ),
        float_format=required_text(output, "float_format"),
        missing_value=required_text(output, "missing_value"),
        create_diagnostic_figure=required_bool(
            diagnostics, "create_diagnostic_figure"
        ),
        show_diagnostic_figure=required_bool(
            diagnostics, "show_diagnostic_figure"
        ),
        copy_namelist_to_output=required_bool(
            output, "copy_namelist_to_output"
        ),
        namelist_copy_filename=required_text(
            output, "namelist_copy_filename"
        ),
    )

    return RunConfiguration(
        source_path=source_path,
        paths=paths,
        controls=controls,
        model=model_options,
        output=output_options,
    )


def _read_csv(path: Path) -> pd.DataFrame:
    """Read a CSV robustly when it may contain a UTF-8 BOM."""

    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")
    frame = pd.read_csv(path, encoding="utf-8-sig")
    frame.columns = [_clean_text(column) for column in frame.columns]
    return frame


def _find_column(
    frame: pd.DataFrame,
    aliases: tuple[str, ...],
    *,
    required: bool = True,
) -> str | None:
    """Find the first matching column name from a list of accepted aliases."""

    normalized = {_clean_text(column).lower(): column for column in frame.columns}
    for alias in aliases:
        match = normalized.get(_clean_text(alias).lower())
        if match is not None:
            return match
    if required:
        raise ValueError(
            "Missing required column. Expected one of: " + ", ".join(aliases)
        )
    return None


def _numeric_series(
    frame: pd.DataFrame,
    aliases: tuple[str, ...],
    *,
    required: bool = True,
    fill_value: float = np.nan,
) -> pd.Series:
    """Load one numeric input column, optionally creating a missing-value series."""

    column = _find_column(frame, aliases, required=required)
    if column is None:
        return pd.Series(fill_value, index=frame.index, dtype=float)
    return pd.to_numeric(frame[column], errors="coerce").astype(float)


def load_meteorology(path: Path, controls: ControlParameters) -> pd.DataFrame:
    """Read meteorological data into a standardized, named-column table."""

    source = _read_csv(path)
    meteorology = pd.DataFrame(index=source.index)
    meteorology["day"] = _numeric_series(source, ("Day", "DOY", "Julian Day"))
    meteorology["hour"] = _numeric_series(source, ("Hour", "Local Hour"))
    meteorology["temperature_c"] = _numeric_series(
        source,
        ("AirTem(degreeC)", "AirTemp(C)", "Temperature(C)", "Tair(C)"),
    )
    meteorology["temperature_k"] = meteorology["temperature_c"] + 273.15
    meteorology["ppfd"] = _numeric_series(
        source,
        ("PPFD(umol/m2/s)", "PPFD", "PAR(umol/m2/s)"),
    )
    meteorology["lai"] = _numeric_series(source, ("LAI", "Leaf Area Index"))
    meteorology["pressure_pa"] = _numeric_series(
        source,
        ("AtmPres(Pa)", "Pressure(Pa)", "PRES"),
        required=controls.humidity_mode == "QV",
    )
    meteorology["wind_m_s"] = _numeric_series(
        source,
        ("WSD(m/s)", "Wind(m/s)", "WindSpeed(m/s)", "WIND"),
    )

    if controls.humidity_mode == "RH":
        meteorology["humidity_input"] = _numeric_series(
            source,
            ("RH(%)", "RH", "RelativeHumidity(%)", "Humidity"),
        )
    else:
        qv_aliases = (
            "QV(kg/kg)",
            "QV",
            "WaterVaporMixingRatio(kg/kg)",
            "SpecificHumidity(kg/kg)",
        )
        qv_column = _find_column(source, qv_aliases, required=False)
        if qv_column is not None:
            meteorology["humidity_input"] = pd.to_numeric(
                source[qv_column], errors="coerce"
            ).astype(float)
        else:
            # The original format reused the fourth column for either RH or QV
            # without necessarily changing its header.  Keep that case usable,
            # but issue an explicit unit warning rather than silently treating
            # percentages as kg kg-1.
            fallback_column = _find_column(source, ("RH(%)", "RH"))
            warnings.warn(
                "RH_QV=0 but no QV-specific column was found. Values in "
                f"'{fallback_column}' will be interpreted as kg kg-1. Rename "
                "the column to QV(kg/kg) to document the units clearly.",
                RuntimeWarning,
                stacklevel=2,
            )
            meteorology["humidity_input"] = pd.to_numeric(
                source[fallback_column], errors="coerce"
            ).astype(float)

    meteorology["isoprene_observed_mg_m2_h"] = _numeric_series(
        source,
        ("Isop(mg/m2/h)", "Isoprene(mg/m2/h)", "Isop_obs"),
        required=False,
    )
    meteorology["swc10_m3_m3"] = _numeric_series(
        source,
        ("SWC10(m3/m3)", "SWC10", "SoilWaterContent10cm"),
        required=False,
    )
    meteorology["kc"] = _numeric_series(
        source,
        ("Kc", "KC"),
        required=False,
    )
    meteorology["kc_7d"] = _numeric_series(
        source,
        ("Kc_7d", "KC_7d", "Kc7d"),
        required=False,
    )
    if meteorology["kc_7d"].isna().all() and not meteorology["kc"].isna().all():
        meteorology["kc_7d"] = meteorology["kc"]

    if meteorology.empty:
        raise ValueError(f"Meteorology file contains no records: {path}")
    if meteorology[["day", "hour"]].isna().any().any():
        raise ValueError("Day and Hour must be present for every meteorological record")
    if not meteorology["hour"].between(0.0, 24.0, inclusive="left").all():
        raise ValueError("Hour values must be in the interval [0, 24)")

    return meteorology


def load_species_parameters(path: Path) -> SpeciesParameters:
    """Read species names, emission factors, and light-dependent fractions."""

    frame = _read_csv(path)
    category_column = _find_column(frame, ("Categories", "Species", "Compound"))
    emission_factor_column = _find_column(
        frame,
        (
            "Emission Factor (nmol/m-2/s-1)",
            "Emission Factor",
            "EF",
        ),
    )
    ldf_column = _find_column(frame, ("LDF", "Light Dependent Fraction"))

    names = tuple(_clean_text(value) for value in frame[category_column])
    emission_factors = pd.to_numeric(
        frame[emission_factor_column], errors="coerce"
    ).to_numpy(dtype=float)
    light_dependent_fraction = pd.to_numeric(
        frame[ldf_column], errors="coerce"
    ).to_numpy(dtype=float)

    if not names:
        raise ValueError("At least one species/activity class is required")
    if len(names) > activity.N_ACTIVITY_CLASSES:
        raise ValueError(
            f"The activity-factor tables support at most "
            f"{activity.N_ACTIVITY_CLASSES} classes; received {len(names)}"
        )
    if not names[0].lower().startswith("isoprene"):
        raise ValueError(
            "The first row of the EF/LDF file must be isoprene because the "
            "current activity-class coefficient ordering assumes index 0."
        )
    if not np.all(np.isfinite(emission_factors)) or np.any(emission_factors < 0.0):
        raise ValueError("Emission factors must be finite and non-negative")
    if not np.all(np.isfinite(light_dependent_fraction)) or np.any(
        (light_dependent_fraction < 0.0) | (light_dependent_fraction > 1.0)
    ):
        raise ValueError("LDF values must be finite and between 0 and 1")

    return SpeciesParameters(names, emission_factors, light_dependent_fraction)


def load_pft_parameters(path: Path) -> PFTParameters:
    """Read and validate the six canopy-type fractions (percent)."""

    frame = _read_csv(path)
    pft_column = _find_column(frame, ("PFT", "Canopy Type", "Vegetation Type"))
    fraction_column = _find_column(frame, ("Fraction", "Fraction(%)", "Percent"))

    names = tuple(_clean_text(value) for value in frame[pft_column])
    fractions = pd.to_numeric(frame[fraction_column], errors="coerce").to_numpy(
        dtype=float
    )

    if len(fractions) != canopy.N_CANOPY_TYPES:
        raise ValueError(
            f"Exactly {canopy.N_CANOPY_TYPES} PFT fractions are required in the "
            f"MEGCAN order; received {len(fractions)}"
        )
    if not np.all(np.isfinite(fractions)) or np.any(fractions < 0.0):
        raise ValueError("PFT fractions must be finite and non-negative")
    if fractions.sum() <= 0.0:
        raise ValueError("At least one PFT fraction must be greater than zero")

    return PFTParameters(names, fractions)


def load_inputs(
    paths: ModelPaths,
    controls: ControlParameters,
) -> InputBundle:
    """Load and validate all model inputs defined by the namelist."""

    meteorology = load_meteorology(paths.meteorology, controls)
    species = load_species_parameters(paths.emission_factors)
    pfts = load_pft_parameters(paths.pft_fractions)
    return InputBundle(meteorology, controls, species, pfts)


def _daily_extreme_by_record(
    days: np.ndarray,
    values: np.ndarray,
    operation: Literal["min", "max"],
) -> np.ndarray:
    """Compute a NaN-safe daily extreme and map it to every record."""

    output = np.full(values.shape, np.nan, dtype=float)
    for day in np.unique(days):
        mask = days == day
        finite = values[mask][np.isfinite(values[mask])]
        if finite.size == 0:
            value = np.nan
        else:
            value = np.min(finite) if operation == "min" else np.max(finite)
        output[mask] = value
    return output


def precompute_meteorology(
    meteorology: pd.DataFrame,
    options: ModelOptions,
) -> PrecomputedMeteorology:
    """Precompute daily statistics once per day rather than once per record."""

    days = meteorology["day"].to_numpy(dtype=float)
    temperature_k = meteorology["temperature_k"].to_numpy(dtype=float)
    wind = meteorology["wind_m_s"].to_numpy(dtype=float)

    grouped = meteorology.groupby("day", sort=False)
    daily_mean_ppfd = grouped["ppfd"].transform("mean").to_numpy(dtype=float)
    daily_mean_temperature = grouped["temperature_k"].transform("mean").to_numpy(
        dtype=float
    )
    previous_10_day_temperature = time_utils.map_previous_daily_mean(
        days,
        temperature_k,
        window_days=10,
    )

    daily_max_temperature = _daily_extreme_by_record(days, temperature_k, "max")
    daily_min_temperature = _daily_extreme_by_record(days, temperature_k, "min")
    daily_max_wind = _daily_extreme_by_record(days, wind, "max")

    return PrecomputedMeteorology(
        daily_mean_ppfd,
        daily_mean_temperature,
        previous_10_day_temperature,
        daily_max_temperature,
        daily_min_temperature,
        daily_max_wind,
    )


def _layer_weights(number_of_layers: int) -> np.ndarray:
    """Return vertical quadrature weights used to integrate activity factors."""

    if number_of_layers == 5:
        return _LAYER_WEIGHTS_5.copy()
    return np.full(number_of_layers, 1.0 / number_of_layers, dtype=float)


def _record_is_valid(
    meteorology: pd.DataFrame,
    record_index: int,
    controls: ControlParameters,
) -> bool:
    """Check whether variables required for canopy calculations are finite."""

    row = meteorology.iloc[record_index]
    required_values = (
        row["day"],
        row["hour"],
        row["ppfd"],
        row["temperature_k"],
        row["humidity_input"],
        row["lai"],
        row["wind_m_s"],
    )
    if controls.humidity_mode == "QV":
        required_values += (row["pressure_pa"],)
    return bool(np.all(np.isfinite(required_values)))


def _lai_pair(
    lai: np.ndarray,
    record_index: int,
) -> tuple[float, float]:
    """Return the previous and current LAI for one record.

    The current record supplies current LAI. If the immediately preceding LAI
    is unavailable, current LAI is used for both values so that one missing
    record does not invalidate the following valid record.
    """

    current_lai = float(lai[record_index])
    if record_index == 0 or not np.isfinite(lai[record_index - 1]):
        previous_lai = current_lai
    else:
        previous_lai = float(lai[record_index - 1])
    return previous_lai, current_lai


def _water_vapor_pressure_pa(
    humidity_input: float,
    temperature_k: float,
    pressure_pa: float,
    controls: ControlParameters,
) -> float:
    """Convert the selected humidity input to canopy vapor pressure (Pa)."""

    if controls.humidity_mode == "RH":
        return canopy.relative_humidity_to_vapor_pressure(
            humidity_input,
            temperature_k,
        )
    return canopy.mixing_ratio_to_vapor_pressure(humidity_input, pressure_pa)


def _calculate_canopy_state(
    *,
    day: float,
    hour: float,
    temperature_k: float,
    ppfd: float,
    lai: float,
    wind_m_s: float,
    humidity_input: float,
    pressure_pa: float,
    controls: ControlParameters,
    pfts: PFTParameters,
    options: ModelOptions,
    layer_positions: np.ndarray,
) -> CanopyState:
    """Calculate the PFT-weighted canopy light and leaf-temperature state."""

    n_layers = options.number_of_layers
    if lai <= 0.0:
        # Emissions are exactly zero after multiplying by LAI.  A finite canopy
        # state keeps intermediate activity-factor calculations well behaved.
        return CanopyState(
            np.full(n_layers, temperature_k, dtype=float),
            np.full(n_layers, temperature_k, dtype=float),
            np.zeros(n_layers, dtype=float),
            np.zeros(n_layers, dtype=float),
            np.zeros(n_layers, dtype=float),
        )

    solar_elevation_deg = time_utils.solar_elevation_angle(
        day,
        controls.latitude_deg,
        hour,
    )
    sin_solar_elevation = np.sin(solar_elevation_deg / 57.29578)
    solar_w_m2 = ppfd / options.solar_to_ppfd
    maximum_solar = (
        sin_solar_elevation
        * options.solar_constant_w_m2
        * time_utils.solar_eccentricity_factor(day)
    )
    (
        diffuse_visible,
        beam_visible,
        diffuse_nir,
        beam_nir,
    ) = time_utils.partition_solar_radiation(solar_w_m2, maximum_solar)

    vapor_pressure_pa = _water_vapor_pressure_pa(
        humidity_input,
        temperature_k,
        pressure_pa,
        controls,
    )

    pft_weights = pfts.fractions_percent / pfts.fractions_percent.sum()
    active_pfts = np.flatnonzero(pft_weights > 0.0)

    sun_temperature = np.zeros(n_layers, dtype=float)
    shade_temperature = np.zeros(n_layers, dtype=float)
    sun_ppfd = np.zeros(n_layers, dtype=float)
    shade_ppfd = np.zeros(n_layers, dtype=float)
    sun_fraction = np.zeros(n_layers, dtype=float)

    # A major performance improvement: canopy radiation/energy balance is
    # calculated only for PFTs with non-zero fractions.  The original code ran
    # all six PFTs even when five had zero weight.
    for pft_index in active_pfts:
        radiation = canopy.canopy_radiation(
            layer_positions,
            n_layers,
            lai,
            sin_solar_elevation,
            beam_visible,
            diffuse_visible,
            beam_nir,
            diffuse_nir,
            int(pft_index),
        )
        temperature_lapse_rate = canopy.canopy_temperature_lapse_rate(
            int(pft_index),
            solar_w_m2,
        )
        energy = canopy.canopy_energy_balance(
            temperature_lapse_rate,
            n_layers,
            layer_positions,
            int(pft_index),
            temperature_k,
            wind_m_s,
            radiation.sun_ppfd,
            radiation.shade_ppfd,
            radiation.sun_visible,
            radiation.shade_visible,
            radiation.sun_nir,
            radiation.shade_nir,
            vapor_pressure_pa,
        )

        weight = pft_weights[pft_index]
        sun_temperature += weight * energy.sun_leaf_temperature_k
        shade_temperature += weight * energy.shade_leaf_temperature_k
        sun_ppfd += weight * radiation.sun_ppfd
        shade_ppfd += weight * radiation.shade_ppfd
        sun_fraction += weight * radiation.sun_fraction

    return CanopyState(
        sun_temperature,
        shade_temperature,
        sun_ppfd,
        shade_ppfd,
        sun_fraction,
    )


def _stress_factors(
    species_index: int,
    *,
    daily_max_temperature_k: float,
    daily_min_temperature_k: float,
    daily_max_wind_m_s: float,
    controls: ControlParameters,
) -> tuple[float, float, float]:
    """Return high-temperature, low-temperature, and high-wind factors."""

    gamma_high_temperature = (
        activity.gamma_ht(species_index, daily_max_temperature_k)
        if controls.high_temperature_response
        else 1.0
    )
    gamma_low_temperature = (
        activity.gamma_lt(species_index, daily_min_temperature_k)
        if controls.low_temperature_response
        else 1.0
    )
    gamma_high_wind = (
        activity.gamma_hw(species_index, daily_max_wind_m_s)
        if controls.high_wind_response
        else 1.0
    )
    return gamma_high_temperature, gamma_low_temperature, gamma_high_wind


def _calculate_record_emissions(
    *,
    current_lai: float,
    previous_lai: float,
    canopy_state: CanopyState,
    daily_mean_ppfd: float,
    daily_mean_temperature_k: float,
    previous_10_day_temperature_k: float,
    daily_max_temperature_k: float,
    daily_min_temperature_k: float,
    daily_max_wind_m_s: float,
    kc_7d: float,
    controls: ControlParameters,
    species: SpeciesParameters,
    options: ModelOptions,
    layer_weights: np.ndarray,
) -> np.ndarray:
    """Calculate all species emissions for one valid model record."""

    n_species = len(species.names)
    emissions = np.full(n_species, np.nan, dtype=float)

    if current_lai <= 0.0:
        emissions.fill(0.0)
        return emissions

    gamma_leaf_age = activity.gamma_a(
        previous_lai,
        current_lai,
        daily_mean_temperature_k,
    )
    gamma_canopy_depth = activity.gamma_cd(
        options.number_of_layers,
        current_lai,
    )
    gamma_bidirectional = (
        activity.gamma_laibidir(current_lai)
        if controls.bidirectional_lai_response
        else 1.0
    )

    # CO2 response is controlled independently by GAMCO2_YN and applies only
    # to isoprene.
    gamma_co2_isoprene = (
        activity.gamma_co2(options.co2_ppm) if controls.co2_response else 1.0
    )

    light_sun = np.array(
        [
            activity.gamp(value, daily_mean_ppfd)
            for value in canopy_state.sun_ppfd
        ],
        dtype=float,
    )
    light_shade = np.array(
        [
            activity.gamp(value, daily_mean_ppfd)
            for value in canopy_state.shade_ppfd
        ],
        dtype=float,
    )

    for species_index in range(n_species):
        light_dependent_fraction = species.light_dependent_fraction[species_index]

        gamma_air_quality = (
            activity.gamma_aq(species_index, options.air_quality_index)
            if controls.air_quality_response
            else 1.0
        )
        (
            gamma_high_temperature,
            gamma_low_temperature,
            gamma_high_wind,
        ) = _stress_factors(
            species_index,
            daily_max_temperature_k=daily_max_temperature_k,
            daily_min_temperature_k=daily_min_temperature_k,
            daily_max_wind_m_s=daily_max_wind_m_s,
            controls=controls,
        )

        if species_index == 0 and controls.soil_moisture_response:
            gamma_soil_moisture = activity.gamma_sm_kc(
                kc_7d,
                options.kc_max,
                options.kc_min,
            )
        else:
            gamma_soil_moisture = 1.0

        gamma_temperature_light = 0.0
        for layer_index in range(options.number_of_layers):
            sun_fraction = canopy_state.sun_fraction[layer_index]
            shade_fraction = 1.0 - sun_fraction

            light_dependent_activity = gamma_canopy_depth[layer_index] * (
                activity.gamtld(
                    canopy_state.sun_leaf_temperature_k[layer_index],
                    daily_mean_temperature_k,
                    previous_10_day_temperature_k,
                    species_index,
                )
                * light_sun[layer_index]
                * sun_fraction
                + activity.gamtld(
                    canopy_state.shade_leaf_temperature_k[layer_index],
                    daily_mean_temperature_k,
                    previous_10_day_temperature_k,
                    species_index,
                )
                * light_shade[layer_index]
                * shade_fraction
            )
            light_independent_activity = (
                activity.gamtli(
                    canopy_state.sun_leaf_temperature_k[layer_index],
                    species_index,
                )
                * sun_fraction
                + activity.gamtli(
                    canopy_state.shade_leaf_temperature_k[layer_index],
                    species_index,
                )
                * shade_fraction
            )
            gamma_temperature_light += layer_weights[layer_index] * (
                light_dependent_activity * light_dependent_fraction
                + light_independent_activity
                * (1.0 - light_dependent_fraction)
            )

        common_activity = (
            current_lai
            * gamma_temperature_light
            * gamma_leaf_age[species_index]
            * gamma_high_wind
            * gamma_air_quality
            * gamma_high_temperature
            * gamma_low_temperature
            * gamma_soil_moisture
        )

        if species_index == 0:
            common_activity *= gamma_co2_isoprene
        elif species_index == 12:
            # The bidirectional LAI factor is applied only to the combined
            # acetaldehyde/ethanol activity class, matching the source code.
            common_activity *= gamma_bidirectional

        emissions[species_index] = (
            common_activity
            * species.emission_factors_nmol_m2_s[species_index]
        )

    return emissions


def simulate(bundle: InputBundle, options: ModelOptions) -> SimulationResult:
    """Run the site-scale MEGAN model and return output tables in memory."""

    meteorology = bundle.meteorology
    controls = bundle.controls
    species = bundle.species
    pfts = bundle.pfts
    precomputed = precompute_meteorology(meteorology, options)

    n_records = len(meteorology)
    n_species = len(species.names)
    emissions = np.full((n_records, n_species), np.nan, dtype=float)

    layer_positions = canopy.gaussian_layer_positions(options.number_of_layers)
    layer_weights = _layer_weights(options.number_of_layers)

    days = meteorology["day"].to_numpy(dtype=float)
    hours = meteorology["hour"].to_numpy(dtype=float)
    temperature_k = meteorology["temperature_k"].to_numpy(dtype=float)
    ppfd = meteorology["ppfd"].to_numpy(dtype=float)
    lai = meteorology["lai"].to_numpy(dtype=float)
    pressure = meteorology["pressure_pa"].to_numpy(dtype=float)
    wind = meteorology["wind_m_s"].to_numpy(dtype=float)
    humidity = meteorology["humidity_input"].to_numpy(dtype=float)
    kc_7d = meteorology["kc_7d"].to_numpy(dtype=float)

    for record_index in range(n_records):
        if not _record_is_valid(
            meteorology,
            record_index,
            controls,
        ):
            continue

        previous_lai, current_lai = _lai_pair(lai, record_index)

        if not np.isfinite(current_lai) or not np.isfinite(previous_lai):
            continue

        canopy_state = _calculate_canopy_state(
            day=days[record_index],
            hour=hours[record_index],
            temperature_k=temperature_k[record_index],
            ppfd=ppfd[record_index],
            lai=current_lai,
            wind_m_s=wind[record_index],
            humidity_input=humidity[record_index],
            pressure_pa=pressure[record_index],
            controls=controls,
            pfts=pfts,
            options=options,
            layer_positions=layer_positions,
        )

        emissions[record_index, :] = _calculate_record_emissions(
            current_lai=current_lai,
            previous_lai=previous_lai,
            canopy_state=canopy_state,
            daily_mean_ppfd=precomputed.daily_mean_ppfd[record_index],
            daily_mean_temperature_k=precomputed.daily_mean_temperature_k[
                record_index
            ],
            previous_10_day_temperature_k=(
                precomputed.previous_10_day_temperature_k[record_index]
            ),
            daily_max_temperature_k=precomputed.daily_max_temperature_k[
                record_index
            ],
            daily_min_temperature_k=precomputed.daily_min_temperature_k[
                record_index
            ],
            daily_max_wind_m_s=precomputed.daily_max_wind_m_s[record_index],
            kc_7d=kc_7d[record_index],
            controls=controls,
            species=species,
            options=options,
            layer_weights=layer_weights,
        )

    all_species = pd.DataFrame(emissions, columns=species.names)
    all_species.insert(0, "hour", hours)
    all_species.insert(0, "day", days)

    nmol_s_to_mg_h = (
        1.0e-9
        * options.isoprene_molecular_weight_g_mol
        * 1000.0
        * 3600.0
    )
    modeled_isoprene = emissions[:, 0] * nmol_s_to_mg_h
    isoprene = pd.DataFrame(
        {
            "day": days,
            "hour": hours,
            "isoprene_observed_mg_m2_h": meteorology[
                "isoprene_observed_mg_m2_h"
            ].to_numpy(dtype=float),
            "isoprene_modeled_mg_m2_h": modeled_isoprene,
            "air_temperature_C": meteorology["temperature_c"].to_numpy(
                dtype=float
            ),
            "ppfd_umol_m2_s": ppfd,
            # Corrected from the original output, which accidentally wrote
            # wind speed under the LAI header.
            "lai_m2_m2": lai,
            "swc10_m3_m3": meteorology["swc10_m3_m3"].to_numpy(dtype=float),
        }
    )

    metrics = calculate_isoprene_metrics(
        isoprene,
        options.daytime_start_hour,
        options.daytime_end_hour,
    )
    return SimulationResult(all_species, isoprene, metrics)


def calculate_isoprene_metrics(
    isoprene: pd.DataFrame,
    daytime_start_hour: float,
    daytime_end_hour: float,
) -> dict[str, float | int | str | None]:
    """Calculate regression and error metrics for daytime isoprene records."""

    daytime = isoprene["hour"].between(
        daytime_start_hour,
        daytime_end_hour,
        inclusive="both",
    )
    observed = isoprene.loc[daytime, "isoprene_observed_mg_m2_h"].to_numpy(
        dtype=float
    )
    modeled = isoprene.loc[daytime, "isoprene_modeled_mg_m2_h"].to_numpy(
        dtype=float
    )
    observed_available = np.isfinite(observed)
    modeled_available = np.isfinite(modeled)
    valid = observed_available & modeled_available

    metrics: dict[str, float | int | str | None] = {
        "daytime_start_hour": float(daytime_start_hour),
        "daytime_end_hour": float(daytime_end_hour),
        # Counts make missing-observation handling explicit in the JSON output.
        # "n" remains the number of paired finite observations used by all
        # error statistics and the regression, preserving the previous API.
        "n_daytime_records": int(daytime.sum()),
        "n_observed_available": int(observed_available.sum()),
        "n_modeled_available": int(modeled_available.sum()),
        "n_missing_observed": int((~observed_available).sum()),
        "n": int(valid.sum()),
        "slope": None,
        "intercept": None,
        "r": None,
        "r_squared": None,
        "p_value": None,
        "rmse_mg_m2_h": None,
        "mae_mg_m2_h": None,
        "mean_bias_mg_m2_h": None,
    }

    observed = observed[valid]
    modeled = modeled[valid]

    if observed.size == 0:
        return metrics

    residual = modeled - observed
    metrics["rmse_mg_m2_h"] = float(np.sqrt(np.mean(residual**2)))
    metrics["mae_mg_m2_h"] = float(np.mean(np.abs(residual)))
    metrics["mean_bias_mg_m2_h"] = float(np.mean(residual))

    if observed.size >= 2 and np.nanstd(observed) > 0.0:
        from scipy import stats

        regression = stats.linregress(observed, modeled)
        metrics["slope"] = float(regression.slope)
        metrics["intercept"] = float(regression.intercept)
        metrics["r"] = float(regression.rvalue)
        metrics["r_squared"] = float(regression.rvalue**2)
        metrics["p_value"] = float(regression.pvalue)

    return metrics


def _all_species_output_table(result: SimulationResult) -> pd.DataFrame:
    """Return the all-species table with a unit in every column header."""

    table = result.all_species.copy()
    if len(table.columns) < 2:
        raise ValueError("The all-species result must contain day and hour columns")

    species_columns = list(table.columns[2:])
    table.columns = [
        "day [day of year]",
        "hour [local decimal hour]",
        *[
            f"{species_name} [nmol m-2 s-1]"
            for species_name in species_columns
        ],
    ]
    return table


def _isoprene_output_table(result: SimulationResult) -> pd.DataFrame:
    """Return the isoprene diagnostic table with explicit physical units."""

    column_names = {
        "day": "day [day of year]",
        "hour": "hour [local decimal hour]",
        "isoprene_observed_mg_m2_h": (
            "isoprene observed [mg m-2 h-1]"
        ),
        "isoprene_modeled_mg_m2_h": "isoprene modeled [mg m-2 h-1]",
        "air_temperature_C": "air temperature [deg C]",
        "ppfd_umol_m2_s": "PPFD [umol m-2 s-1]",
        "lai_m2_m2": "LAI [m2 m-2]",
        "swc10_m3_m3": "soil water content at 10 cm [m3 m-3]",
    }
    missing_columns = [
        column for column in column_names if column not in result.isoprene.columns
    ]
    if missing_columns:
        raise ValueError(
            "The isoprene result is missing expected columns: "
            + ", ".join(missing_columns)
        )
    return result.isoprene.rename(columns=column_names)


def write_outputs(
    result: SimulationResult,
    output_directory: Path,
    *,
    output_options: OutputOptions | None = None,
    namelist_path: Path | None = None,
) -> dict[str, Path]:
    """Write tables, metrics, and an optional copy of the run namelist.

    CSV headers are generated at write time so the in-memory tables retain
    compact programmatic names while every file column exposes its unit.
    """

    settings = output_options or OutputOptions()
    output_directory.mkdir(parents=True, exist_ok=True)
    all_species_path = output_directory / settings.all_species_filename
    isoprene_path = output_directory / settings.isoprene_filename
    metrics_path = output_directory / settings.metrics_filename
    for path in (all_species_path, isoprene_path, metrics_path):
        path.parent.mkdir(parents=True, exist_ok=True)

    _all_species_output_table(result).to_csv(
        all_species_path,
        index=False,
        float_format=settings.float_format,
        na_rep=settings.missing_value,
    )
    _isoprene_output_table(result).to_csv(
        isoprene_path,
        index=False,
        float_format=settings.float_format,
        na_rep=settings.missing_value,
    )
    with metrics_path.open("w", encoding="utf-8") as file_handle:
        json.dump(result.metrics, file_handle, indent=2, ensure_ascii=False)
        file_handle.write("\n")

    output_paths = {
        "all_species": all_species_path,
        "isoprene": isoprene_path,
        "metrics": metrics_path,
    }

    if settings.copy_namelist_to_output and namelist_path is not None:
        namelist_copy_path = output_directory / settings.namelist_copy_filename
        namelist_copy_path.parent.mkdir(parents=True, exist_ok=True)
        source = namelist_path.expanduser().resolve()
        destination = namelist_copy_path.resolve()
        if source != destination:
            shutil.copy2(source, destination)
        output_paths["namelist"] = namelist_copy_path

    return output_paths


def save_diagnostic_figure(
    result: SimulationResult,
    output_path: Path,
    *,
    daytime_start_hour: float = 9.0,
    daytime_end_hour: float = 17.0,
    show: bool = False,
) -> Path:
    """Save a time-series and observed-vs-modeled isoprene diagnostic figure."""

    import matplotlib.pyplot as plt

    isoprene = result.isoprene
    daytime = isoprene["hour"].between(
        daytime_start_hour,
        daytime_end_hour,
        inclusive="both",
    )
    subset = isoprene.loc[daytime]
    observed = subset["isoprene_observed_mg_m2_h"].to_numpy(dtype=float)
    modeled = subset["isoprene_modeled_mg_m2_h"].to_numpy(dtype=float)

    # Convert hour to a fraction of a day; the original code added raw hours to
    # DOY, which distorted the time coordinate.
    time_coordinate = (
        subset["day"].to_numpy(dtype=float)
        + subset["hour"].to_numpy(dtype=float) / 24.0
    )

    figure, axes = plt.subplots(1, 2, figsize=(12, 5.5), constrained_layout=True)
    axes[0].plot(time_coordinate, observed, "o", markersize=3, label="Observation")
    axes[0].plot(time_coordinate, modeled, "o", markersize=3, label="Model")
    axes[0].set_xlabel("Day of year")
    axes[0].set_ylabel(r"Isoprene flux (mg m$^{-2}$ h$^{-1}$)")
    axes[0].legend()

    valid = np.isfinite(observed) & np.isfinite(modeled)
    if np.any(valid):
        axes[1].plot(observed[valid], modeled[valid], "o", markersize=4)
        maximum = float(np.nanmax(np.concatenate([observed[valid], modeled[valid]])))
        upper = max(1.0, maximum * 1.05)
        axes[1].plot([0.0, upper], [0.0, upper], "-")
        axes[1].set_xlim(0.0, upper)
        axes[1].set_ylim(0.0, upper)
    axes[1].set_xlabel(r"Observed isoprene (mg m$^{-2}$ h$^{-1}$)")
    axes[1].set_ylabel(r"Modeled isoprene (mg m$^{-2}$ h$^{-1}$)")

    slope = result.metrics.get("slope")
    intercept = result.metrics.get("intercept")
    r_squared = result.metrics.get("r_squared")
    if slope is not None and intercept is not None and r_squared is not None:
        axes[1].text(
            0.04,
            0.96,
            f"y = {slope:.2f}x + {intercept:.2f}\n$R^2$ = {r_squared:.2f}",
            transform=axes[1].transAxes,
            verticalalignment="top",
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=200)
    if show:
        plt.show()
    plt.close(figure)
    return output_path
