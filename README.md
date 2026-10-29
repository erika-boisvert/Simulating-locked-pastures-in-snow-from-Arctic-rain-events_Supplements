# Simulating-locked-pastures-in-snow-from-Arctic-rain-events_Supplements
Supplements for Boisvert-Vigneault et al., paper on using SVS2-Crocus to simulate locked pastures from rain-on-snow events in the Arctic.

## Meteorological datasets
### ERA5
[Link to download](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=download) \
[Link to documentation](https://confluence.ecmwf.int/spaces/CKB/pages/76414402/ERA5+data+documentation)

### CaSRv3.2
[Link to download](https://pavics.ouranos.ca/datasets.html#c) \
[Link to documentation](https://hpfx.collab.science.gc.ca/~scar700/rcas-casr/dataset_specifics.html#available_variables)

### Environment and Climate Change Canada Weather Stations
*Station IDs*:
  - Sachs Harbour: 1794 (from 1980 to 2015), 53326 (from 2015 onwards) and 53328 (from 2015 onwards with weather observations)
  - Thomsen River: 10598
  - Nangmagvik Lake: 27118
\
[Link to download](https://climate.weather.gc.ca/historical_data/search_historic_data_e.html) \
The R-script "download_bulk_ECCC-LWS.R" allows to download larger timeseries at once and look up stations by name or location.

## SVS2-Crocus
[Link to GitHub official release](https://github.com/VVionnet/MESH_SVS/releases/tag/v1.0.0) \
