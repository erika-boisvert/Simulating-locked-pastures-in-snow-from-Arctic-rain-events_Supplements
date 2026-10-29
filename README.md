# Simulating-locked-pastures-in-snow-from-Arctic-rain-events_Supplements
Supplements for Boisvert-Vigneault et al., paper on using SVS2-Crocus to simulate locked pastures from rain-on-snow events in the Arctic.

## Abstract
Rain-on-snow (ROS) events occurring in the Arctic can form dense ice layers, preventing ungulates such as muskoxen and caribou from foraging throughout the winter. Such snow conditions are known as locked pastures as the vegetation is trapped beneath the dense ice layer. These locked pastures can have severe consequences for Arctic wildlife populations and Indigenous communities that depend on them. This study (1) assessed how well two meteorological reanalyses, ERA5 and the Canadian Surface Reanalysis (CaSR) v3.2, represent the meteorological conditions driving ROS events, and (2) evaluated whether the multi-physics snowpack model Crocus can reproduce the resulting locked-pasture stratigraphy. The reanalyses were evaluated against three local weather stations on Banks Island in the Canadian Arctic, for near-surface air temperature, wind speed, relative humidity, and precipitation. Four methods to determine the precipitation phase were also compared against weather observations. CaSRv3.2 paired with an atmospheric-profile-based partitioning method reproduced the timing and phase of a ROS event in 2003 most accurately and was used to drive SVS2-Crocus in an ensemble composed of four water holding capacity equations, three snow-compaction schemes, three falling snow density representations, and three wind packing variations. The ensemble reproduced a locked pasture consistent with reported snow conditions in 79 of 108 members of the ensemble for the 2003–2004 winter. Sensitivity experiments that removed either the rainfall or the accompanying warm temperatures showed that above-freezing temperatures driving melt-freeze events were the dominant driver of snowpack densification, although rainfall increased the probability of a locked pasture. A 35-year (1985–2020) simulation validated against remote sensing ROS detections identified accurately four other winters that were subject to a ROS event in late fall that resulted in a dense basal snow layer. These results show that physically based snowpack models can capture and help explain the formation of ROS-driven locked pastures once meteorological reanalysis biases are identified.


## Datasets
### Meteorological datasets
#### ERA5
[Link to download](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=download) \
[Link to documentation](https://confluence.ecmwf.int/spaces/CKB/pages/76414402/ERA5+data+documentation)

#### CaSRv3.2
[Link to download](https://pavics.ouranos.ca/datasets.html#c) \
[Link to documentation](https://hpfx.collab.science.gc.ca/~scar700/rcas-casr/dataset_specifics.html#available_variables)

#### Environment and Climate Change Canada Weather Stations
*Station IDs*:
- Sachs Harbour: 1794 (from 1980 to 2015), 53326 (from 2015 onwards) and 53328 (from 2015 onwards with weather observations)
- Thomsen River: 10598
- Nangmagvik Lake: 27118 \
\
[Link to download](https://climate.weather.gc.ca/historical_data/search_historic_data_e.html) \
The R-script "download_bulk_ECCC-LWS.R" allows to download larger timeseries at once and look up stations by name or location.

#### Adjusted Daily Rainfall and Snowfall Dataset for Canada
[Link to download](https://open.canada.ca/data/en/dataset/d8616c52-a812-44ad-8754-7bcc0d8de305) \
[Link to documentation](https://www.tandfonline.com/doi/full/10.1080/07055900.2017.1342163)

### SVS2-Crocus Snowpack Model
[Link to GitHub official release](https://github.com/VVionnet/MESH_SVS/releases/tag/v1.0.0) \
[Link to documentation](https://mesh-model.atlassian.net/wiki/spaces/USER/overview?homepageId=6389894)

## Structure of the repository
**SVS2-Crocus_figures/fig07-10_normalized_ensemble_350.py** \
*Python code to recreate figure 07 and figure 10(a) and (b) from SVS2-Crocus netCDF outputs.* \
\
**SVS2-Crocus_figures/fig09_Simulation-4panels.py** \
*Python code to recreate figure 09 from SVS2-Crocus netCDF outputs.* \
\
**SVS2-Crocus_figures/fig11_timeseries_ensemble_350_barplot.py** \
*Python code to recreate figure 11 from SVS2-Crocus netCDF outputs.* \
\
**SVS2-Crocus_params/Crocus_ensemble.sh** \
*Bash script used to run SVS2-Crocus as an ensemble over varying parameters.* \
\
**SVS2-Crocus_params/MESH_input_run_options.ini** \
*Initiation script to specify the input format of the soil and meteorological driving data for SVS2-Crocus.* \
\
**SVS2-Crocus_params/MESH_input_soil_levels.txt** \
*Text file identifying the soil layering for SVS2-Crocus.* \
\
**SVS2-Crocus_params/MESH_parameters_orig.ini** \
*Initiation script that specifies all options within the snowpack model. This is the file that Crocus_ensemble.sh iterates through to change the snowpack schemes for the ensemble.* \
\
**SVS2-Crocus_params/generate_nc_output.py** \
*Python script to create a netCDF file out of SVS2-Crocus standard generated outputs. It is used at the end of Crocus_ensemble.sh* \
\
**code_for_figures.ipynb** \
*Jupyter Notebook containing the codes to recreate most of the figures from the paper, expect for SVS2-Crocus output figures. Most input files required to run this notebook are issued from "meteorological_data_processing.ipynb"; simulation specific files are outputs of SVS2-Crocus_figures data processing.* \
\
**meteorological_data_processing.ipynb** \
*Jupyter Notebook containing the codes to process and analyse the meteorological data, as done in the paper. The output files are directly used to create the figures in "code_for_figures.ipynb".* \
\
**download_bulk_ECCC-LWS.R** \
*R script used to download the Environment and Climate Change Canada weather station data in bulk.* \
