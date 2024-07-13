### Correlation Coefficients

- `calc_CC` : Contains function to calculates correlation coefficient matrix for a 1D time-series
 

- `cc_APX`: Example CC calculation for a timeseries of random noise using
  the `calculate_CC` function defined in `calc_CC.py`


- `cc_earth`: Computes CC matrix and mean cc for Earth timeseries. Saves plot to `figs`.
 

- `cc_lunar`: Computes CC matrix and mean cc for Lunar timeseries, as well as a small submatrix 
              between user-defined frequency bands. See Fig 9 of paper. Saves plot to `figs`.


- `cc_mars`: Computes CC matrix and mean cc for Martian timeseries. Saves plot to `figs`.


- `cc_ps`: Computes CC matrix and mean cc for Parameter Space timeseries. Saves data to user defined
  path currently set to `data/CC/T=10_PROPERCUTOFF`. 