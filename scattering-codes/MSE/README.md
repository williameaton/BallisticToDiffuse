### Multi scale entropy: 

All scripts are in relevant directories. Each directory is listed below with a description
of the scripts. 


#### matlab 
- Contains some old matlab scripts written and provided by Dr Surya Pacchai. 
- Ultimately not used in this study (since python MSE available) but maybe
  of use to someone


#### mmse 
**TODO: move MMSE up so its in same 
directory as planetary etc and make subsection here for it
Currently in planetraey**


#### parameter_space 
- `PS_MSE`: Computes MSE for parameter space timeseries for the 200 seconds afer
  maximum energy arrival. 

#### planetary 
*Note that there is a separate MSE file for each data type since sometimes specific glitches etc need to be removed*

- `earth_mse`: Computes MSE and integrated MSE of Earth samples. Saves MSE to file. Also computes power spectra.


- `lunar_mse`: Computes MSE and integrated MSE of Earth samples. Saves MSE to file.  Also computes power spectra.


- `mars_mse_S0173`: Computes MSE and integrated MSE of Mars event S0173. Saves MSE to file.
                    Slices out glitches in timeseries and saves the deglitched timeseries for MMSE analysis.


- `mars_mse_S0235`: Same as `mars_mse_S0173` but for S0235. Note that these scripts could definitely be
  combined into a single function -- I was clearly very lazy when this was coded! 


- `plot_freq_planet_comparison`: Reproduces Fig A13 (supplementary material) plot of moon/mars MSE at different 
  bandpasses. Need to run the lunar and mars MSE scripts first. 


- `ps_example`: Example MSE calculation for PS data



#### synthetic_mse_apx 

- `checking_MSE`: 


- `example_apx_MSE_seis`: 


- `plot_data`
