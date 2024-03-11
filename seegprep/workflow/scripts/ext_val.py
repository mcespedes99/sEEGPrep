import pyedflib
import pandas as pd
import numpy as np
import re
from pandas.api.types import is_numeric_dtype
import sys
import logging
import json


# All channels from the channels.tsv should be on the EDF
# getPhysicalDimensions from EDF should match the units from channels.tsv
# filters should match between files (low, high and notch) edf, channels.tsv (not the json bc names of the filters don't have a defined format). (COMPLETE)
# srate should match between files edf, channels and json (COMPLETE)
def check_channels(edf_obj, tsv_obj, json_obj):
    # Get info from channels.tsv (only SEEG, ECOG or EEG)
    channels_tsv = (
        tsv_obj[tsv_obj["type"] == "SEEG"]["name"].to_list()
        + tsv_obj[tsv_obj["type"] == "EEG"]["name"].to_list()
        + tsv_obj[tsv_obj["type"] == "ECOG"]["name"].to_list()
    )
    # Get channels from edf
    channels_edf = edf_obj.getSignalLabels()
    units_edf = [header["dimension"] for header in edf_obj.getSignalHeaders()]
    srate_edf = edf_obj.getSampleFrequencies() / edf_obj.datarecord_duration
    # Check all from channels are in edf
    labels_warnings = []
    units_warnings = []
    srate_warnings = []
    filters_warnings = []
    for idx_tsv, channel in enumerate(channels_tsv):
        # Check labels
        if channel not in channels_edf:
            labels_warnings.append(
                f"Warning: Channel {channel} is present on the channels.tsv file but not on the corresponding EDF file."
            )
        else:
            idx_edf = channels_edf.index(channel)
            # Check units
            if units_edf[idx_edf] != tsv_obj.loc[idx_tsv, "units"].replace("μ", "u"):
                units_warnings.append(
                    f"Warning: Channel {channel} has different units on the channels.tsv file ({tsv_obj.loc[idx_tsv,'units']}) compared to the EDF file ({units_edf[idx_edf]})."
                )
            # Check srate
            if srate_edf[idx_edf] != tsv_obj.loc[idx_tsv, "sampling_frequency"]:
                srate_warnings.append(
                    f"Warning: Channel {channel} has different sampling rate on the channels.tsv file ({tsv_obj.loc[idx_tsv,'sampling_frequency']}) compared to the EDF file ({srate_edf[idx_edf]})."
                )
            # Check filters
            list_filt = [
                ("Low Pass", "low_cutoff", "LP"),
                ("High Pass", "high_cutoff", "HP"),
                ("Notch", "notch", "N"),
            ]
            for filt_name, filt_tsv, filt_edf in list_filt:
                # Check for LP in the edf file
                pattern = f"{filt_edf}:(\d+\.*\d*)Hz"
                string = edf_obj.getPrefilter(idx_edf)
                if not (
                    np.isnan(tsv_obj.loc[idx_tsv, filt_tsv])
                    or (filt_edf != "LP" and tsv_obj.loc[idx_tsv, filt_tsv] == 0)
                ):
                    # Case 1: If it's set in the channels.tsv and not on the edf file
                    if not re.search(pattern, string):
                        filters_warnings.append(
                            f"Warning: {filt_name} filter defined in the channels.tsv for the channel {channel} but not in the EDF file."
                        )
                    # Case 2: Different values
                    elif (
                        re.search(pattern, string).group(1)
                        != tsv_obj.loc[idx_tsv, filt_tsv]
                    ):
                        filters_warnings.append(
                            f"Warning: The {filt_name} filter for the channel {channel} has different cutoff frequency on the channels.tsv file ({tsv_obj.loc[idx_tsv, filt_tsv]}) compared to the one defined on the EDF file ({re.search(pattern, string).group(1)})."
                        )
                # Case 3: Defined on edf but not on tsv
                elif (np.isnan(tsv_obj.loc[idx_tsv, filt_tsv])) and (
                    re.search(pattern, string)
                ):
                    filters_warnings.append(
                        f"Warning: {filt_name} filter defined in the EDF for the channel {channel} but not in the channels.tsv file."
                    )
    # Check srate against json file using edf file
    if not (srate_edf == json_obj["SamplingFrequency"]).all():
        channels_edf = np.array(channels_edf)
        chns_problem = channels_edf[~(srate_edf == json_obj["SamplingFrequency"])]
        if len(chns_problem) > 1:
            srate_warnings.append(
                f"Warning: Channels {chns_problem} have different sampling rate on the EDF file ({ srate_edf[~(srate_edf == json_obj['SamplingFrequency'])]}) compared to the json file ({json_obj['SamplingFrequency']})."
            )
        else:
            srate_warnings.append(
                f"Warning: Channel {chns_problem} has different sampling rate on the EDF file ({ srate_edf[~(srate_edf == json_obj['SamplingFrequency'])][0]}) compared to the json file ({json_obj['SamplingFrequency']})."
            )
    # Check srate against json file using tsv file
    srate_tsv = tsv_obj["sampling_frequency"].to_numpy()
    if not (srate_tsv == json_obj["SamplingFrequency"]).all():
        chns_problem = np.array(channels_tsv)[
            ~(srate_tsv == json_obj["SamplingFrequency"])
        ]
        if len(chns_problem) > 1:
            srate_warnings.append(
                f"Warning: Channels {chns_problem} have different sampling rate on the channels.tsv file ({ srate_tsv[~(srate_tsv == json_obj['SamplingFrequency'])]}) compared to the json file ({json_obj['SamplingFrequency']})."
            )
        else:
            srate_warnings.append(
                f"Warning: Channel {chns_problem} has different sampling rate on the channels.tsv file ({ srate_tsv[~(srate_tsv == json_obj['SamplingFrequency'])][0]}) compared to the json file ({json_obj['SamplingFrequency']})."
            )
    # Check results
    if len(labels_warnings) == 0:
        labels_warnings.append(
            f"Info: Test PASSED. All channels present on the channels.tsv file were found in the corresponding EDF file."
        )
    if len(units_warnings) == 0:
        units_warnings.append(
            f"Info: Test PASSED. Units in the channels.tsv file match the units in the EDF file."
        )
    if len(srate_warnings) == 0:
        srate_warnings.append(
            f"Info: Test PASSED. Sampling frequency in the files channels.tsv, EDF and json file match."
        )
    if len(filters_warnings) == 0:
        filters_warnings.append(
            f"Info: Test PASSED. Filters defined in the channels.tsv match with the defined filters in the EDF file."
        )
    return labels_warnings + units_warnings + filters_warnings + srate_warnings


# x,y,z should be numbers on the electrodes.tsv
def check_format_electrodes(electrodes_tsv):
    pos_warnings = []
    if not (
        is_numeric_dtype(electrodes_tsv["x"])
        and is_numeric_dtype(electrodes_tsv["y"])
        and is_numeric_dtype(electrodes_tsv["z"])
    ):
        pos_warnings.append(
            "Warning: Electrodes positions in the electrodes.tsv file are not numeric, as required by the BIDS standard."
        )
    if len(pos_warnings) == 0:
        pos_warnings.append(
            "Info: Test PASSED. All electrodes positions in the electrodes.tsv file are numeric, as required by the BIDS standard."
        )
    return pos_warnings


# All contacts from channels in the edf and channels should be present on the electrodes.tsv
def check_electrodes(
    tsv_obj, elec_obj, json_obj, modality
):  # Reference should be bipolar or unipolar
    ref_dict = {"ieeg": "iEEGReference", "eeg": "EEGReference"}
    miss_electrodes = []
    assumption = False
    if json_obj[ref_dict[modality]] == "bipolar":
        miss_electrodes.append(
            "Info: Using bipolar reference based on the parameter 'iEEGReference' from the json file."
        )
        # Extract labels
        channels_channelstsv = tsv_obj["name"].tolist()
        electrodes = elec_obj["name"].tolist()
        # Only checking channels.tsv as we already checked that all the channels from this file are in the EDF based on the function check_channels
        # The EDF file might have some other channels with no position
        for chn in channels_channelstsv:
            # Extract unipolar labels
            pattern = r"([a-zA-Z]+(\d+-)?)(\d+)-+(\d+)"
            try:
                matches = re.search(pattern, chn).groups()
                unipolar_chn1 = matches[0] + matches[2]
                unipolar_chn2 = matches[0] + matches[3]
                if unipolar_chn1 not in electrodes:
                    miss_electrodes.append(
                        f"Warning: Electrode {unipolar_chn1} found in channels.tsv file but not in electrodes.tsv file."
                    )
                if unipolar_chn2 not in electrodes:
                    miss_electrodes.append(
                        f"Warning: Electrode {unipolar_chn2} found in channels.tsv file but not in electrodes.tsv file."
                    )
            except:
                raise TypeError("Channels do not seem to have a bipolar format")
    else:
        assumption = True
        miss_electrodes.append(
            "Info: Assuming unipolar reference as the parameter 'iEEGReference' is not set to 'bipolar' in the json file."
        )
        # Extract labels
        channels_channelstsv = tsv_obj["name"].tolist()
        electrodes = elec_obj["name"].tolist()
        for chn in channels_channelstsv:
            if chn not in electrodes:
                miss_electrodes.append(
                    f"Warning: Electrode {chn} found in channels.tsv file but not in electrodes.tsv file."
                )
    if (len(miss_electrodes) == 1 and assumption) or len(miss_electrodes) == 0:
        miss_electrodes.append(
            "Info: Test PASSED. All the contacts forming the channels in the channels.tsv and EDF files are present in the electrodes.tsv."
        )
    return miss_electrodes


# number of channels per category should match the number of channels in edf and channels.tsv
def check_number_channels(edf_obj, tsv_obj, json_obj):
    # Check number of channels per category for channels tsv file
    valid_types = ["EEG", "EOG", "ECG", "EMG", "ECOG", "SEEG"]
    chn_numbers_warnings = []
    for type_mod in valid_types:
        key_json = f"{type_mod}ChannelCount"
        length_chns = tsv_obj[tsv_obj.type == type_mod].shape[0]
        if key_json in json_obj:
            if length_chns != json_obj[key_json]:
                chn_numbers_warnings.append(
                    f"Warning: json file indicates that {json_obj[key_json]} {type_mod} channels should be present but {length_chns} {type_mod} channels were found in the channels.tsv file."
                )
    # Check only 'EEG', 'SEEG' and 'ECOG' againts EDF file
    brain_types = ["EEG", "ECOG", "SEEG"]
    total_brain_chns = 0
    for type_mod in valid_types:
        key_json = f"{type_mod}ChannelCount"
        if key_json in json_obj:
            total_brain_chns += json_obj[key_json]

    total_chns_edf = len(edf_obj.getSignalLabels())
    if total_chns_edf < total_brain_chns:
        chn_numbers_warnings.append(
            f"Warning: json file indicates that {total_brain_chns} EEG+SEEG+ECOG channels should be present but only {total_chns_edf} channels were found in the EDF file."
        )

    if len(chn_numbers_warnings) == 0:
        chn_numbers_warnings.append(
            "Info: Test PASSED. Number of channels specified in the json file match the number of channels in the channels.tsv and EDF files."
        )
    return chn_numbers_warnings


# Checks that the duration of the file in the json file matches the one on the edf file
def check_duration(edf_obj, json_obj, tolerace=0.1):
    import decimal

    # Tolerace in seconds
    # Get duration from edf
    edf_dur = np.round(edf_obj.file_duration, 1)
    # Get duration from json_obj
    json_dur = json_obj["RecordingDuration"]
    # Check if they match
    file_duration_warnings = []
    if (
        not np.round(json_dur - tolerace, 1)
        <= edf_dur
        <= np.round(json_dur + tolerace, 1)
    ):
        file_duration_warnings.append(
            f"Warning: Length specified in the edf file ({edf_dur} seconds) doesn't seem to match the length indicated in the JSON file ({json_dur})."
        )
    # Find number of decimals in json_dur
    d = decimal.Decimal(f"{json_dur}")
    n_dec = abs(int(d.as_tuple().exponent))
    # Try to check if the json file has the length specified in other units
    possible_cases = [(60, "minutes"), (3600, "hours")]
    for t, case in possible_cases:
        tmp_edf_dur = np.round(edf_dur / t, n_dec)
        if json_dur - (tolerace / t) <= tmp_edf_dur <= json_dur + (tolerace / t):
            file_duration_warnings.append(
                f"Warning: It seems that the JSON file is currently indicating the file length in {case}, but the standard defines RecordingDuration in seconds."
            )

    return file_duration_warnings


def main():
    edf = snakemake.input.edf
    channels_tsv = snakemake.input.channels_tsv
    tsv_elec = snakemake.input.tsv_elec
    json_file = snakemake.input.json_file
    modality = snakemake.params.modality
    out_txt = snakemake.output.out_txt
    LOG_FILENAME = snakemake.log[0]
    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG)
    try:
        # Read EDF
        edf_obj = pyedflib.EdfReader(edf)
        # Read channels tsv
        chn_tsv_obj = pd.read_csv(channels_tsv, sep="\t")
        # Read electrodes tsv
        elec_obj = pd.read_csv(tsv_elec, sep="\t")
        # Read json
        with open(json_file, "r") as f:
            json_obj = json.load(f)
        # Start running validation
        bids_val_results = []
        bids_val_results += check_channels(edf_obj, chn_tsv_obj, json_obj)
        bids_val_results += check_format_electrodes(elec_obj)
        bids_val_results += check_electrodes(chn_tsv_obj, elec_obj, json_obj, modality)
        bids_val_results += check_number_channels(edf_obj, chn_tsv_obj, json_obj)
        bids_val_results += check_duration(edf_obj, json_obj)
        # Save into txt obj
        with open(out_txt, "w") as outfile:
            outfile.writelines([line + "\n" for line in bids_val_results])

    except:
        logging.exception("Got exception on main handler")
        raise


if __name__ == "__main__":
    main()
