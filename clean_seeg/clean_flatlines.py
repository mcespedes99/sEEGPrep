"""
Clean flatlines from sEEG file.
This code is was recovered from the repository
https://github.com/moeinrazavi/EEG-ASR-Python
which was originally inspired in the Matlab implementation found in 
https://github.com/sccn/clean_rawdata
Author of EEG-ASR-Python: Moein Razavi
Author of this repository: Mauricio Cespedes Tenorio
"""
#!/usr/bin/env python
# coding: utf-8

import numpy as np


def clean_flatlines(signal,SRate,MaxFlatlineDuration=5,MaxAllowedJitter=20):
    
    # Remove (near-) flat-lined channels.
    # signal = clean_flatlines(signal,MaxFlatlineDuration,MaxAllowedJitter)
    #
    # This is an automated artifact rejection function which ensures that 
    # the data contains no flat-lined channels.
    #
    # In:
    #   signal (ndarray) : continuous data set, assumed to be appropriately 
    #                      high-passed (e.g. >0.5Hz or
    #                      with a 0.5Hz - 2.0Hz transition band). 
    #                      Structure: (n_epochs x) n_channels x samples
    #                      n_epochs is optional
    #
    #   MaxFlatlineDuration : Maximum tolerated flatline duration. In seconds. If a channel has a longer
    #                         flatline than this, it will be considered abnormal. Default: 5
    #
    #   MaxAllowedJitter : Maximum tolerated jitter during flatlines. As a multiple of epsilon.
    #                      Default: 20
    #
    # Out:
    #   signal : data set with flat channels removed
    #
    # Examples:
    #   % use with defaults
    #   eeg = clean_flatlines(eeg);
    ## TODO: include epochs management!

    if signal.ndim == 3:
        n_epochs, n_channels, n_time = np.shape(signal)
        raise Exception('Only 2D arrays at supported for now.')
    elif signal.ndim == 3:
        n_channels, n_time = np.shape(signal)
        n_epochs = 1
    else:
        raise Exception('The signal must be at least 2D (n_channels x samples).')

    # flag channels
    removed_channels = np.array([False for i in range(len(n_channels))])
    

    for c in range(n_channels):
        
        zero_intervals = np.reshape(np.where(np.diff([False] + list(abs(np.diff(signal[c,:]))<(MaxAllowedJitter*np.finfo(float).eps))+[False])), (2,-1)).T
        
        if (len(zero_intervals) > 0):
            if (np.max(zero_intervals[:,1] - zero_intervals[:,0]) > MaxFlatlineDuration*SRate):
                removed_channels[c] = True
    new_channels_inds = np.where(~removed_channels)
    # remove them
    if all(removed_channels):
        print('Warning: all channels have a flat-line portion; not removing anything.')
    elif any(removed_channels):
        print('Now removing flat-line channels...')
        
        signal = signal[new_channels_inds]
    return signal, new_channels_inds[0]
