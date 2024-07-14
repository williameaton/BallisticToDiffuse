# Imports
import numpy as np

def gen_mmse_slices(trace, slice_length, overlap, normalise=True):
    # Splits up a timeseries into slices for moving MSE calculation
    # Function is currently used and script written so that it may be easily called in a shell script but can be edited
    # slightly to work as a single script, so long as variables are defined in-script.

    # Calculating number of slices
    data_len = len(trace)
    no_slices = int(np.floor((data_len-overlap)/(slice_length - overlap)))

    ## Initialise slices array:
    slices_array = np.zeros((no_slices, slice_length))

    # Normalise the trace:
    if normalise:
        trace_norm = trace/np.amax(np.abs(trace))
    else:
        trace_norm = trace

    # Slice up the array:
    for sl in range(no_slices):
        slices_array[sl, :] = trace_norm[ (slice_length-overlap)*sl:((slice_length-overlap)*sl)+slice_length]

    print(f"generated {no_slices} slices")
    return slices_array, no_slices

