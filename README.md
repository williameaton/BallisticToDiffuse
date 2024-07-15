### About: 

This repository contains all codes used in the paper:

[**Seismic scattering regimes from multiscale entropy and frequency correlations**](https://doi.org/10.1093/gji/ggae098)

- Data for synthetic seismograms are hosted at ... due to filesize issues.
- Planetary data can be created by running the download scripts in the ```data_acquisition_processing``` directory. 

Authors: Will Eaton, Claudia Haindl, Tarje Nissen-Meyer \
Contact: weaton@princeton.edu 

Last updated: July 11 2024



### Data access: 

Accompanying these scripts is ~ 150 Gb of data. Evidently this can not be hosted in this GitHub repository and so is hosted at the following DOI: 

doi.org/10.5281/zenodo.12745350

This tarball can be untarred using 

```bash
  $ tar -xvf BallisticTODiffuseData.tar.gz 
```

The relative paths used in the scripts rely on the `Data` directory being within the main parent directory. Once extracted, your main directory should have the following structure: 

├── data \
├── README.md\
├── plotting\
├── scattering-codes\
├── simulation_metadata\
└── wetools.py
