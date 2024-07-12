#!/bin/bash

#pip install obspy
# USER DEFINE VARIABLES:
# First define the variables for downloading the mars data:
#data_fname='S0173a'
#stime='2019-05-23T02:15:00.000'
#etime='2019-05-23T03:15:00.000'

data_fname='S0235b'
stime='2019-07-26T12:05:00.000'
etime='2019-07-26T13:20:00.000'

# Download data and process using python function:
python3 mars_download.py $data_fname $stime $etime
