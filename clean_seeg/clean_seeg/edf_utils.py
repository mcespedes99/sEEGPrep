import re
import numpy as np 
import pyedflib
from .val_utils import (
    get_chn_labels
)


def extract_signal(channels_tsv, edf_path):
    # Retrive channels from channels.tsv file that are present in the EDF file
    chn_labels, discarded_labels = get_chn_labels(channels_tsv, elec_edf)
    with pyedflib.EdfReader(edf_path) as edf_in:
        # Extract labels from EDF file
        elec_edf = edf_in.getSignalLabels()
        N = edf_in.getNSamples()[0]
        signal = []
        for chan in chn_labels:
            id_ch = elec_edf.index(chan)
            chn_sig = edf_in.readSignal(id_ch)
            signal.append(chn_sig)
        signal = np.vstack(signal)
    return chn_labels, discarded_labels, signal

def adjust_n_samples(n_samples, f):
    remainder = n_samples % f.getSampleFrequencies()[0]
    if remainder != 0:
        n_samples += f.getSampleFrequencies()[0] - remainder
    return n_samples

def find_timestamps(event_start, f, srate, n_samples):
    t = np.arange(0, f.getNSamples()[0]) / srate
    if event_start is not None:
        if isinstance(event_start, int):
            time_stamps_init = t[event_start]
        else:
            id = [
                value[0]
                for value in enumerate(f.readAnnotations()[2])
                if re.match(event_start, value[1], re.IGNORECASE)
            ]
            onset_list = f.readAnnotations()[0]
            time_stamps_init = onset_list[id]
        if not hasattr(time_stamps_init, '__iter__'):
            time_stamps_init = [time_stamps_init]
    else:
        time_stamps_init = t[::n_samples]
    return time_stamps_init