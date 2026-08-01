# MEGANv3 Site-Scale Model in Python

This repository is a standalone site-scale implementation of the Model of
Emissions of Gases and Aerosols from Nature (MEGAN). It estimates emissions for
19 BVOC/CO activity classes from meteorological drivers, canopy composition,
emission factors, and environmental activity factors.

The model has one supported calculation workflow. All run settings are stored
in `MEGAN_namelist.ini`; normal users do not edit `main_program.py`.

## Run

From the project directory:

```bash
python main_program.py
```

Use another namelist without modifying the launcher:

```bash
python main_program.py --namelist experiments/site_b.ini
```

Relative input and output paths are resolved from the directory containing the
selected namelist.

## Project structure

```text
MEGAN3_structured/
├── main_program.py              # only executable model launcher
├── MEGAN_namelist.ini           # all run settings and file paths
├── inputs/                      # model input data only
│   ├── Met_HourlyData_2012_moflux_Kc.csv
│   ├── EF_LDF.csv
│   └── PFT_Fraction.csv
├── src/                         # model source code
│   ├── __init__.py
│   ├── megan_model.py
│   ├── MEGVEA.py
│   ├── MEGCAN.py
│   └── TIMEFUNC.py
├── outputs/                     # generated model products
├── tests/                       # automated tests
├── PROJECT_OVERVIEW_CN.md
├── VALIDATION.md
├── requirements.txt
├── environment.yml
└── LICENSE
```
`main_program.py` is intentionally kept at the project root as the user-facing
launcher. All reusable model code is imported from the `src` package, while all
CSV input datasets are isolated in `inputs`. Normal users therefore edit only
`MEGAN_namelist.ini` and files under `inputs/`.

## Namelist

### `[FILES]`

```ini
[FILES]
meteorology = inputs/Met_HourlyData_2012_moflux_Kc.csv
emission_factors = inputs/EF_LDF.csv
pft_fractions = inputs/PFT_Fraction.csv
output_directory = outputs
```

### `[MODEL]`

```ini
[MODEL]
number_of_layers = 5
solar_constant_w_m2 = 1361.5
solar_to_ppfd = 2.1
air_quality_index = 40.0
co2_ppm = 403.0
kc_min = 0.0
kc_max = 0.82
isoprene_molecular_weight_g_mol = 68.12
```

There is no mode selector. The model always uses the maintained calculation
logic, including current-record LAI, NaN-safe daily extrema, and the independent
`GAMCO2_YN` switch.

### `[SITE]`

```ini
[SITE]
RH_QV = 1
Latitude = 45.4
WT = 0.196
```

- `RH_QV = 1` or `RH`: relative humidity in percent.
- `RH_QV = 0` or `QV`: water-vapor mixing ratio in kg kg⁻¹; atmospheric
  pressure is then required.
- `Latitude`: decimal degrees; north is positive.
- `WT`: soil wilting point used by the soil-moisture response function.

### `[ACTIVITY_SWITCHES]`

```ini
[ACTIVITY_SWITCHES]
GAMBD_YN = 0
GAMAQ_YN = 0
GAMHT_YN = 0
GAMLT_YN = 0
GAMHW_YN = 0
GAMSM_YN = 0
GAMCO2_YN = 0
```

Boolean values accept `0/1`, `true/false`, `yes/no`, or `on/off`.

### `[DIAGNOSTICS]` and `[OUTPUT]`

These sections define the daytime comparison window, figure behavior, output
filenames, numeric precision, missing-value text, and whether the exact
namelist is copied into the output directory.

## Meteorological input

The bundled file is `inputs/1.Met_HourlyData_2012_moflux_Kc.csv`. Variables are read
by header name, not by fixed column position.

| Variable | Example header | Unit | Requirement |
|---|---|---:|---|
| Day of year | `Day` | 1–366 | Required |
| Local hour | `Hour` | decimal hour | Required |
| Air temperature | `AirTem(degreeC)` | °C | Required |
| Relative humidity | `RH(%)` | % | Required when `RH_QV=1` |
| Water-vapor mixing ratio | `QV(kg/kg)` | kg kg⁻¹ | Required when `RH_QV=0` |
| PPFD | `PPFD(umol/m2/s)` | µmol m⁻² s⁻¹ | Required |
| LAI | `LAI` | m² m⁻² | Required |
| Atmospheric pressure | `AtmPres(Pa)` | Pa | Required when `RH_QV=0` |
| Wind speed | `WSD(m/s)` | m s⁻¹ | Required |
| Isoprene observation | `Isop(mg/m2/h)` | mg m⁻² h⁻¹ | Optional |
| Soil water at 10 cm | `SWC10(m3/m3)` | m³ m⁻³ | Optional |
| ET/PET ratio | `Kc` | dimensionless | Optional |
| Seven-day ET/PET ratio | `Kc_7d` | dimensionless | Needed for drought response |

### Missing isoprene observations

`Isop(mg/m2/h)` is diagnostic only. Individual values may be blank or `NaN`,
and the entire column may be omitted. Missing observations do not affect model
emissions. They remain missing in the diagnostic CSV and are excluded from
paired statistics and the observation-versus-model scatter plot.

## Outputs

The default output directory contains:

- `Moflux_2012_simulation_all.csv`
- `Moflux_2012_simulation_isoprene.csv`
- `Moflux_2012_metrics.json`
- `Moflux_2012_isoprene_diagnostics.png`
- `MEGAN_namelist_used.ini`

The all-species CSV reports every emission class in `nmol m-2 s-1`. The
isoprene diagnostic CSV reports observed and modeled isoprene in
`mg m-2 h-1`. Every CSV header contains its physical unit or coordinate
meaning.

The metrics JSON records the number of daytime records, available observations,
available simulations, paired records, regression statistics, RMSE, MAE, and
mean bias. It contains no calculation-mode field.

## Tests

Run:

```bash
python -m unittest discover -s tests -v
```

The tests cover the namelist, unit-bearing output headers, optional isoprene
observations, maintained LAI handling, NaN-safe daily extrema, independent CO₂
response, RH/QV humidity input, and all six canopy types under day and night
conditions.

## Reference

Guenther, A. B., Jiang, X., Heald, C. L., Sakulyanontvittaya, T., Duhl, T.,
Emmons, L. K., & Wang, X. (2012). The Model of Emissions of Gases and Aerosols
from Nature version 2.1 (MEGAN2.1): an extended and updated framework for
modeling biogenic emissions. *Geoscientific Model Development*, 5, 1471–1492.

Wang, H., Lu, X., Seco, R., Stavrakou, T., Karl, T., Jiang, X., et al. (2022).
Modeling isoprene emission response to drought and heatwaves within MEGAN using
evapotranspiration data and by coupling with the Community Land Model.
*Journal of Advances in Modeling Earth Systems*, 14, e2022MS003174.
