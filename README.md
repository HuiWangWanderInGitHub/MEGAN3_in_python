# MEGANv3_in_python
## _Discription_
The Model of Emissions of Gases and Aerosols from Nature (MEGAN)(Guenther et al., 2012) is a widely used model framework for estimating BVOC fluxes from terrestrial ecosystem. The MEGANv3 in python is a stand alone version MEGAN designed for the site scale flux estimation. The code is written in python, and it is easy to use with the basic knowledge of Python. The main purpose of this version MEGAN is for model tuning and development.

## _Model Usage_

### _Prerequisite Python packages_
The code package includes the code files as well as input files. The following external Python packages are required for running the model:
* [numpy](https://numpy.org/)
* [pandas](https://pandas.pydata.org/)
* [matplotlib](https://matplotlib.org/)
* [scipy](https://scipy.org/)

Before running the code, please ensure the above packages have been added to your enviroment.
### _Code Files_
There are four code files in the package:
* **_main_program.py_**: The file contains the main program. 
* **_MEGVEA\.py_**: The file contains the functions or algorithms for calculating the impact of different environmental factors (e.g., light and temperature) for BVOC emission.
* **_MEGCAN\.py_**: The file contains the functions for the canopy model in MEGAN.
* **_TIMEFUNC.py_**: The file contains the other trivia functions (e.g., the function for calculating average temperature).

Users can run the model by execute the **_main_program.py_** file. For instance, under the command line environment, users can run the model by typing:

`python3   ./main_program.py`

Users who uses the platform like [Anaconda](https://www.anaconda.com/) can just click the Running button to execute the main program to run the MEGAN. The input files and ouput files are set up in the first few codes in "main_program.py" as:

`Inputfile = "./1.Met_HourlyData_2012_moflux_Kc.csv"`

`ParaInputFile = "./2.Parameter.csv"`

`EmisInputFile = "./3.EF_LDF.csv"`

`PFTInputFile = "./4.PFT_Fraction.csv"`

`Output_name_all = 'Moflux_2012_simulaiton.csv'`

`Output_name_isop = 'Moflux_2012_simulation_isop.csv'`

Users can change the name of the input and output files by modifying the above lines in "main_program.py". The description of the files are included in this document.

The scientific algorithms for describing the response of BVOC to different environmental factors are presented in **_MEGCAN\.py_**. Users can modify or improve the algorithms based on their purposes.

Canopy model related functions are in **_MEGCAN\.py_**, and other small functions for the model are in **_TIMEFUNC.py_**. Users in general will keep these functions same unless they are interested to improve the canopy model in MEGAN.
### _Input Files_
The input files in the package contains model running parameters and driven fields. All input files are in _**CSV**_ format and can be easily editted with Excel. There are four input files:
* 1.Met_HourlyData_2012_moflux_Kc.csv
* 2.Parameter.csv
* 3.EF_LDF.csv
* 4.PFT_Fraction.csv

#### _1.Met_HourlyData_2012_moflux_Kc.csv_
"1.Met_HourlyData_2012_moflux_Kc.csv" stores the input fields includes:
| Field | Description | Unit | Comment|
| :---        |    :----:   |          :---: | :---:  |
| Day   | Julian Day  | [1-366] |---|
| Hour   | Hour  |  [0-23]  |---|
| AirTem   | Air Temperature  |  °C  |---|
| PPFD   | Photosynthetic Photon Flux Density |  umol/m2/s  |---|
| LAI  | Leaf Area Index  |  m2/m2  |---|
| RH  | Relative Humidty  |  %  |---|
| AtmPres | Atmospheric Pressure  |  Pa  |---|
| WSD | Wind Speed  |  m/s  |---|
| Isop | Isoprene Observation  |  mg/m2/h  |optional|
| SWC10 | Soil Water Content at 10 cm  |  m3/m3  |optional|
| KC | The ratio of the actual evoportranspiration (ET) to the potential evapotranspiration (PET)  | ---|optional|
| KC_7d | 7-day running averaged ratio of ET to PET  | --- |optional|

When users prepare their own input files, please follow the same format as our example in the package, and missing values can be blank in the csv file. The  the _optional_ field in the "1.Met_HourlyData_2012_moflux_Kc.csv" should be blank if it is not available but with the headers there.
#### _2.Parameter.csv_
"2.Parameter.csv" contains the parameters controlling the stressed algorithms and other model running related parameters. MEGANv3 considers five environmental stresses including:
* Heatwave
* Coldwave
* Air Quality
* High Wind
* Drought

These algorithms are still in a developing stage, so we recommend turning off these algorithms (set the value to 0) when you use MEGANv3 currently. The other parameters in "2.Parameter.csv" are listed as below:

| Field | Description |  Comment|
| :---        |    :----:   |          :---: |
| RH_QV   | 1 for Relative Humidity; 0: water vapor mixing ratio  |---|
| Latitude   | Latitude of the site  |---|
| WT   | Wilting point of the soil   | For the drought algorithm in Guenther et al., 2012|

"_RH\_QV_" is for deciding what kind of the air humidity metric will be used in "1.Met_HourlyData_2012_moflux_Kc.csv". If Water vapor mixing ratio were adopted, users should change "_RH\_QV_" from 1 to 0.

#### _3.EF_LDF.csv and 4.PFT_Fraction.csv_
"3.EF_LDF.csv" saves the Emission Factor (EF, nmol/m2/s) and Light Depedent Fraction (LDF) of the whole ecosystem. "4.PFT_Fraction.csv" stores the fractions of the Plant Functional Types (PFTs). There are 6 types of PFT in the canopy model of MEGANv3:
* Needleleaf Trees
* Tropical Trees
* Temperate Broadleaf Trees
* Shrubs
* Herbaceous
* Crop

Users can modify the emission factors, LDF and fraction of PFTs based on their site environment.
## _Output Files_
The model has two ouput files:

`Output_name_all = 'Moflux_2012_simulaiton.csv'`

`Output_name_isop = 'Moflux_2012_simulation_isop.csv'`

Users can modify the name of the files in "main_program.py". "Output_name_all" is the ouput contains all 19 compounds with the unit of nmol/m2/s in MEGANv3 and "Output_name_isop" is the ouput only for isoprene with the unit of "ug/m2/h"

## _Reference_
Guenther, A. B., Jiang, X., Heald, C. L., Sakulyanontvittaya, T., Duhl, T., Emmons, L. K., & Wang, X. (2012). The Model of Emissions of Gases and Aerosols from Nature version 2.1 (MEGAN2. 1): an extended and updated framework for modeling biogenic emissions. Geoscientific Model Development, 5(6), 1471-1492.
Wang, H., Lu, X., Seco, R., Stavrakou, T., Karl, T., Jiang, X., ... & Guenther, A. B. (2022). Modeling isoprene emission response to drought and heatwaves within MEGAN using evapotranspiration data and by coupling with the community land model. Journal of advances in modeling earth systems, 14(12), e2022MS003174.


