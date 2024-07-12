### Data acquisition & Processing

This directory contains a number of scripts used to download and process data.
Below, I group them into scripts for Earth, Mars, Lunar and synthetic data. 


#### Earth 
- `download_earth`: Downloads data for terrestrial seismograms used in study. 
  Data are filtered, rotated to NEZ and have their instrument response removed. 


- `rotate_earth`: Loads processed Earth data stored in `data/earth/processed` as 
  ZNE and rotates to RTZ, stored in `data/processed_rotated`


#### Lunar 
- `process_lunar_data`: Processes local lunar data stored in `data/lunar` by removing
  instrument response and bandpassing if desired. 

#### Mars
- `download_mars`: Downloads data for martian seismograms used in study. 
  Data are filtered, rotated to NEZ and have their instrument response removed. 


- `rotate_mars`: Rotates processed martian data stored in `data/mars/filtered` rotates
  to RTZ (stored in same directory but has letter 'r' indicating rotation)


- `run_mars_download.sh`: Bash script can be used to download mars data if lines 63-72 of
  `download_mars` are not commented out. 

#### Synthetic 
- `pick_arrival_P`: Reads in synthetic data from the `data/processed/` directory and 
  selects the P wave arrival time based on an amplitude threshold. Pick is saved to 
  the `data/picks/p_arrival` directory. 


- `pick_max_amplitude`: Reads in synthetic data from the `data/processed` directory and 
  selects the time of maximum energy. Pick is saved to the `data/picks/maxE_pick` directory. 


- `process_raw_data`: Decimates/resamples and bandpasses the raw synthetic data (`data/raw`) and
  outputs to `data/processed`


- `slice_by_time`: Processes synthetic data from `data/processed/` and slices it. It slices the 
  time-series starting at the maximum energy pick value `data/picks/maxE_arrival` and ending 
  200 seconds later. The data is outputted to `data/MSE_Sslices_200`. 


#### Other 
- `rotate`: Some functions to rotate Obspy traces/streams. 


