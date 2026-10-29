# Simulating-locked-pastures-in-snow-from-Arctic-rain-events_Supplements
Supplements for Boisvert-Vigneault et al., paper on using SVS2-Crocus to simulate locked pastures from rain-on-snow events in the Arctic.

## Structure of the repository


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
