import pyedflib
import numpy as np
import os, fnmatch
from collections import OrderedDict
from pathlib import Path
from typing import Union
from functools import partial
import traceback
from multiprocessing.pool import Pool
import re
import pandas as pd
import nibabel as nb
import mne

# from EDFlib.edfreader import EDFreader
# from EDFlib.edfwriter import EDFwriter
import sys
import SimpleITK as sitk


# Function to create bipolar channels from given unipolars
def create_bipolar_comb(id, dict_key, channels_dict):
    # Only create bipolar channels with the closest or second closest electrode
    if (
        abs(int(channels_dict[dict_key][id]) - int(channels_dict[dict_key][id + 1]))
        <= 2
    ):
        bipolar_chn = (
            dict_key
            + channels_dict[dict_key][id]
            + "-"
            + channels_dict[dict_key][id + 1]
        )
        chn1_label = dict_key + channels_dict[dict_key][id]
        chn2_label = dict_key + channels_dict[dict_key][id + 1]
        return (bipolar_chn, chn1_label, chn2_label)


def create_bipolar_combi(id, bipolar_list):
    return bipolar_list[id]


# Function to discard
def apply_bipolar_criteria(df_bipolar, bipolar_list, processes):
    non_white_matter_unknown_bool = (
        df_bipolar["region name"].str.contains(
            "White-Matter|Unknown", case=False, regex=True
        )
        == False
    )
    ids = df_bipolar.loc[
        non_white_matter_unknown_bool, "region name"
    ].index.values.tolist()
    # print(ids)
    # Extract bipolar list (original channels + bipolar channel)
    with Pool(processes=processes) as pool:
        filtered_list = pool.map(
            partial(create_bipolar_combi, bipolar_list=bipolar_list), ids
        )
    return filtered_list


# Function to extract position of bipolar channels
def bipolar_info(id, dict_key, channels_dict, elec_pos, df_cols):
    # Only create bipolar channels with the closest or second closest electrode
    if (
        abs(int(channels_dict[dict_key][id]) - int(channels_dict[dict_key][id + 1]))
        <= 2
    ):
        # Bipolar channel label
        bipolar_chn = (
            dict_key
            + channels_dict[dict_key][id]
            + "-"
            + channels_dict[dict_key][id + 1]
        )
        # Channels' labels
        chn1_label = dict_key + channels_dict[dict_key][id]
        chn2_label = dict_key + channels_dict[dict_key][id + 1]
        # Extract positions
        inf_chn1 = elec_pos.loc[elec_pos[df_cols["label"]] == chn1_label]
        inf_chn2 = elec_pos.loc[elec_pos[df_cols["label"]] == chn2_label]
        # print(inf_chn1)
        data = {
            "type": inf_chn1[df_cols["type"]].values[0],
            "group": inf_chn1[df_cols["group"]].values[0],
            "label": bipolar_chn,
            "x": (inf_chn1[df_cols["x"]].values[0] + inf_chn2[df_cols["x"]].values[0])
            / 2,
            "y": (inf_chn1[df_cols["y"]].values[0] + inf_chn2[df_cols["y"]].values[0])
            / 2,
            "z": (inf_chn1[df_cols["z"]].values[0] + inf_chn2[df_cols["z"]].values[0])
            / 2,
        }
        return data


# Function to create a bipolar channel list from
def create_bipolars(electrodes_df, electrodes_edf, processes, df_cols=None):
    # df_cols (dict) = {type_record, label, x, y, z, group}
    # the dict can be in any order. The df will have the structure given
    # by the dict
    if df_cols == None:
        df_cols = {
            "type": "type",
            "label": "label",
            "x": "x",
            "y": "y",
            "z": "z",
            "group": "group",
        }
    channels = {}
    # Try to options of labels
    pattern1 = r"([A-Z0-9]+[-]+)(\d+)$"  # for electrodes like 'LOpS-10' or 'LOpS1-10'
    pattern2 = r"([A-Z]+[-]*)(\d+)$"  # for electrodes like 'LOpS10'
    # Extract channels info
    for electrode in electrodes_df[df_cols["label"]].values:
        if electrode in electrodes_edf:
            match = re.match(pattern1, electrode, re.IGNORECASE)
            if not match:
                match = re.match(pattern2, electrode, re.IGNORECASE)
            if match.group(1) in channels.keys():
                channels[match.group(1)].append(match.group(2))
            else:
                channels[match.group(1)] = [match.group(2)]
    # Create new list
    bipolar_list = []
    bipolar_info_dicts = []
    for key in channels.keys():
        # Extract bipolar list (original channels + bipolar channel)
        with Pool(processes=processes) as pool:
            bipolar_list = bipolar_list + pool.map(
                partial(create_bipolar_comb, dict_key=key, channels_dict=channels),
                list(range(len(channels[key]) - 1)),
            )
        # Extract bipolar channels information
        with Pool(processes=processes) as pool:
            bipolar_info_dicts = bipolar_info_dicts + pool.map(
                partial(
                    bipolar_info,
                    dict_key=key,
                    channels_dict=channels,
                    elec_pos=electrodes_df,
                    df_cols=df_cols,
                ),
                list(range(len(channels[key]) - 1)),
            )
    # Remove None elements
    bipolar_info_dicts = [
        element for element in bipolar_info_dicts if element is not None
    ]
    bipolar_list = [element for element in bipolar_list if element is not None]
    # Convert dict to DataFrame
    bipolar_elec = pd.DataFrame(columns=list(df_cols.keys()))
    data = pd.DataFrame(bipolar_info_dicts)
    bipolar_elec = pd.concat([bipolar_elec, data], ignore_index=True)
    return bipolar_list, bipolar_elec


# Function to extract info from each channel
def extract_channel_data(chn_number, edf_file, chn_list):
    edf_in = pyedflib.EdfReader(edf_file)
    # Get labels from original edf file
    channels_labels = edf_in.getSignalLabels()
    # print('labels')
    # print(channels_labels)
    if type(chn_list[0]) == list or type(chn_list[0]) == tuple:  # Case 1: Bipolar
        # Get indexes of channels
        # print('new iter')
        chn1_id = channels_labels.index(chn_list[chn_number][1])
        chn2_id = channels_labels.index(chn_list[chn_number][2])
        # print(chn_number)
        signal_chn1 = edf_in.readSignal(chn1_id)
        signal_chn2 = edf_in.readSignal(chn2_id)
        chn_data = signal_chn1 - signal_chn2
        # Deallocate space in memory
        edf_in.close()
        del signal_chn1
        del signal_chn2
    elif type(chn_list[0]) == str:  # Case 2: unipolar
        # Get indexes of channel
        # print('new iter')
        chn_id = channels_labels.index(chn_list[chn_number])
        # print(chn_number)
        chn_data = edf_in.readSignal(chn_id)
        # Deallocate space in memory
        edf_in.close()
    return chn_data


# Function to extract headers for bipolar channels
def extract_channel_header(chn_number, original_headers, chn_list, channels_labels):
    if (
        type(chn_list[0]) == list or type(chn_list[0]) == tuple
    ):  # Case 1: Rereferencing to bipolar
        # Get indexes of channels
        chn1_id = channels_labels.index(chn_list[chn_number][1])
        # Update header
        chn_header = original_headers[chn1_id]
        chn_header["label"] = chn_list[chn_number][0]
    elif type(chn_list[0]) == str:  # Case 2: no rereferencing
        # Get indexes of channel
        chn_id = channels_labels.index(chn_list[chn_number])
        # Copy header
        chn_header = original_headers[chn_id]
    return chn_header


# Function to get label map
def get_colors_labels(colortable_file):
    colortable = pd.read_table(colortable_file, index_col="index")
    cols = colortable.columns
    if ("r" in cols) and ("g" in cols) and ("b" in cols):
        colortable[["r", "g", "b"]] = colortable[["r", "g", "b"]].astype("int64")
    else:
        # Generate random colors
        rgb = set()
        while len(rgb) < len(colortable.index):
            r, g, b = (
                np.random.randint(0, 256),
                np.random.randint(0, 256),
                np.random.randint(0, 256),
            )
            if (r, g, b) not in rgb and (r, g, b) != (
                0,
                0,
                0,
            ):  # (0,0,0) reserved for unknown
                rgb.add((r, g, b))
        r, g, b = zip(*rgb)
        # Add to dataframe
        colortable["r"] = r
        colortable["g"] = g
        colortable["b"] = b
    # Append element '0' (unknown) if doesn't exist
    if 0 not in colortable.index:
        unknown_dict = {"name": ["Unknown"], "r": [0], "b": [0], "g": [0]}
        for column in cols:
            if column not in unknown_dict:
                unknown_dict[column] = [""]
        # Concatenate dataframes
        tmp_df = pd.DataFrame(unknown_dict)
        colortable = pd.concat([tmp_df, colortable])
    return colortable


# Function to get label id based on parcellation obj
# TODO: test function in jupyter notebook
def get_electrodes_id(parc, elec_df, df_cols, tfm_list):
    # Load data of parcellations
    data_parc = np.asarray(parc.dataobj)
    # Coordinates in MRI RAS
    mri_ras_mm = elec_df[[df_cols["x"], df_cols["y"], df_cols["z"]]].values
    # print(mri_ras_mm)
    # Apply transforms
    for tfm, inv_bool in tfm_list:
        if type(tfm) == str:
            if tfm.endswith("txt"):
                tfm = readRegMatrix(tfm)
                if inv_bool:
                    tfm = np.linalg.inv(tfm)
                mri_ras_mm = mne.transforms.apply_trans(tfm, mri_ras_mm)
            elif tfm.endswith("nii.gz"):
                # reads the transform and casts the output compaitble format
                transform_image = sitk.ReadImage(tfm)
                transform_image = sitk.Cast(transform_image, sitk.sitkVectorFloat64)
                # load it as a transform
                identity_transform = sitk.Transform(transform_image)
                # Convert points from RAS to LPS
                mri_mni_lps = mri_ras_mm * np.array([-1, -1, 1])
                # Transform
                for point_id in range(mri_mni_lps.shape[0]):
                    mri_mni_lps[point_id, :] = np.array(
                        identity_transform.TransformPoint(mri_mni_lps[point_id, :])
                    )
                # Convert from LPS back to RAS
                mri_ras_mm = mri_mni_lps * np.array([-1, -1, 1])
        else:
            if inv_bool:
                tfm = np.linalg.inv(tfm)
            mri_ras_mm = mne.transforms.apply_trans(tfm, mri_ras_mm)
    # Transform from contrast mri ras to non-contrast MRI ras
    # mri_ras_mm = mne.transforms.apply_trans(non_cont_to_cont_tf, mri_ras_mm)
    # Update position of electrodes in df
    elec_df[df_cols["x"]] = mri_ras_mm[:, 0]
    elec_df[df_cols["y"]] = mri_ras_mm[:, 1]
    elec_df[df_cols["z"]] = mri_ras_mm[:, 2]
    # print(mri_ras_mm)
    # To voxels
    inv_affine = np.linalg.inv(parc.affine)
    # here's where the interpolation should be performed!!
    vox = (mne.transforms.apply_trans(inv_affine, mri_ras_mm)).astype(int)
    # print(vox)
    # Try to get all the indexes
    # For hippunfold, the space is cropped, so it might be the case that one of the vox values is outside of the cropped space.
    # Assign unknown in this case.
    if (
        (vox[:, 0] >= data_parc.shape[0]).any()
        or (vox[:, 1] >= data_parc.shape[1]).any()
        or (vox[:, 2] >= data_parc.shape[2]).any()
    ):
        id = []
        for idx in range(vox.shape[0]):
            # If outside of the cropped space:
            if (
                (vox[idx, 0] >= data_parc.shape[0])
                or (vox[idx, 1] >= data_parc.shape[1])
                or (vox[idx, 2] >= data_parc.shape[2])
            ):
                id.append(0)
            else:
                id.append(data_parc[vox[idx, 0], vox[idx, 1], vox[idx, 2]])
        id = np.array(id)
    # Get the indexes more efficiently
    else:
        id = data_parc[vox[:, 0], vox[:, 1], vox[:, 2]]
    return id, elec_df


# Function to get rgb values for each contact
def get_label_rgb(parc, elec_df, tfm_list, label_map, df_cols):
    # vox, data_parc = ras2vox(parc, elec_df, non_cont_to_cont_tf)
    id, elec_df = get_electrodes_id(parc, elec_df, df_cols, tfm_list)
    print(label_map)
    print(id)
    print(label_map.loc[id, ["name", "r", "g", "b"]])
    vals = label_map.loc[id, ["name", "r", "g", "b"]].to_numpy()
    vals = np.c_[id, vals]
    vals = pd.DataFrame(data=vals, columns=["region ID", "region name", "r", "g", "b"])
    vals = pd.concat([elec_df, vals], axis=1)
    return vals


# Function to read matrix
def readRegMatrix(trsfPath):
    with open(trsfPath) as (f):
        return np.loadtxt(f.readlines())


# Function to create tsv with bipolar channels info
def extract_location(parc_path, chn_info_df, df_cols, tfm_list, colortable_file):
    import os

    # Load labels from LUT file
    labels = get_colors_labels(colortable_file)
    # Load parcellation file
    parc_obj = nb.load(parc_path)
    # Create df
    df = get_label_rgb(parc_obj, chn_info_df, tfm_list, labels, df_cols)
    # if not os.path.exists(out_tsv_name):
    #     df.to_csv(out_tsv_name, sep = '\t')
    return df


# Function to extract useful information from csv file
def get_chn_info(csv_file, electrodes_edf, df_cols=None):  # , conf = 'unipolar'
    df = pd.read_csv(csv_file, sep="\t")
    # df_cols (dict) = {type_record, label, x, y, z, group}
    # the dict can be in any order. The df will have the structure given
    # by the dict
    if df_cols == None:
        df_cols = {
            "type": "type",
            "label": "label",
            "x": "x",
            "y": "y",
            "z": "z",
            "group": "group",
        }
    # print(list(df_cols.values()))
    important_data = df[list(df_cols.values())]
    important_data.reset_index()  # make sure indexes pair with number of rows
    elec_df = pd.DataFrame(columns=list(df_cols.keys()))
    for index in range(len(important_data)):
        if important_data.loc[index, df_cols["label"]] in electrodes_edf:
            tmp_df = pd.DataFrame(
                [important_data.loc[index, list(df_cols.keys())].values],
                columns=list(df_cols.keys()),
            )
            elec_df = pd.concat([elec_df, tmp_df], axis=0)
    del tmp_df
    # reset index
    elec_df = elec_df.reset_index(drop=True)
    # if conf == 'unipolar':
    chn_list = elec_df[df_cols["label"]].values.tolist()
    # elif conf == 'bipolar':
    #     # All the labels must be bipolar!!
    #     labels = df[df_cols['label']].values.tolist()
    #     pattern = r'^(\D+)(\d+)-(\d+)$'
    #     chn_list = []
    #     for label in labels:
    #         groups_list = list(re.search(pattern=pattern, string=label).groups())
    #         first_chn = groups_list[0]+groups_list[1]
    #         second_chn = groups_list[0]+groups_list[2]
    #         full_bip = re.search(pattern=pattern, string=label).group()
    #         chn_list.append([full_bip, first_chn, second_chn])
    return elec_df, chn_list, df_cols


# Function to establish sampling rate of EDF file
def get_srate_params(srate, max_srate=254):  # max_srate=254 as uses int8
    # Convert sampling rate to integer
    srates = [srate, int(srate)]  # try first on srate, if not, on the int
    for rate in srates:
        n_min = int(np.floor(rate / max_srate))
        for i in range(n_min + 1, int(rate / 2)):
            if rate % i == 0:
                duration = 1 / i
                if len(str(duration).split(".")[1]) <= 4:
                    return 1 / i, rate / i
    raise Exception("No appropriate parameters found in get_srate_params")


# Function to create EDF file based in channel list
def create_EDF(
    edf_file,
    out_path,
    processes,
    chn_labels=None,
    signal=None,
    n_removed=None,
    new_srate=None,
):
    import re

    # If signal == None, data is extracted from edf file based on the chn_labels
    try:
        edf_in = pyedflib.EdfReader(edf_file)
        # First import labels
        labels = edf_in.getSignalLabels()
        # Create file:
        if chn_labels == None:
            chn_labels = labels
        edf_out = pyedflib.EdfWriter(
            out_path, len(chn_labels), file_type=pyedflib.FILETYPE_EDFPLUS
        )
        # First set the data from the header of the edf file:
        edf_out.setHeader(edf_in.getHeader())
        headers_orig = edf_in.getSignalHeaders()

        # Set each channel info:
        # if signal == None:
        N = edf_in.getNSamples()[0]
        # else:
        #     N = len(signal[0]) # 0 index chosen randomly
        # Sampling rate:
        if new_srate == None:
            # f.datarecord_duration gives the value is sec and setDatarecordDuration receives it in units
            # of 10 ms. Therefore: setDatarecordDuration = datarecord_duration*10^6 / 10
            # This actually is used to set the sample frequency as the max number that can be written in headers is 254 (uses int8)
            edf_out.setDatarecordDuration(int(edf_in.datarecord_duration * 100000))
        srate = edf_in.getSampleFrequencies()[0] / edf_in.datarecord_duration

        # Annotations
        annot_orig = edf_in.readAnnotations()
        # Fill n_removed if not indicated
        # if n_removed == None:
        #     n_removed = np.zeros(len(annot_orig)+1)

        # Close file
        edf_in.close()
        print(f"len annot {len(annot_orig)}")
        # print(f'len n_rem {len(n_removed)}')
        # Write annotations
        t = np.arange(0, N + 1) / srate
        # The plus 1 is required since the t indicates start points, for the end epoch, it's an end point of
        # the last datarecord!
        # print(t[-10:])
        # pattern_start = r'Epoch #\d starts.'
        # pattern_end = r'Epoch #\d ends.'
        # pattern_id = 0
        for annot_id in np.arange(len(annot_orig)):
            # print(annot_id)
            # print(annot_orig[0][annot_id])
            t_id = np.abs(np.subtract(t, annot_orig[0][annot_id])).argmin()
            t_annot = t[t_id]
            edf_out.writeAnnotation(
                t_annot, annot_orig[1][annot_id], annot_orig[2][annot_id]
            )
            # This was written for the autoreject case, which is currently not working
            # if re.match(pattern_start, annot_orig[2][annot_id]):
            #     edf_out.writeAnnotation(max(0, annot_orig[0][annot_id]-(n_removed[annot_id]*srate)),
            #                             annot_orig[1][annot_id], annot_orig[2][annot_id])
            # elif re.match(pattern_end, annot_orig[2][annot_id]):
            #     edf_out.writeAnnotation(max(0, annot_orig[0][annot_id]-(n_removed[annot_id+1]*srate)),
            #                             annot_orig[1][annot_id], annot_orig[2][annot_id])
            # else:
            #     edf_out.writeAnnotation(max(0, annot_orig[0][annot_id]), annot_orig[1][annot_id], annot_orig[2][annot_id])

        # Extract channel information:
        chn_lists = range(len(chn_labels))
        # print(chn_labels)
        # print(chn_lists)
        print("Channel part")
        if signal == None:  # extracting data from edf file
            # Create bipolar signals:
            # print(chn_lists)
            with Pool(processes=processes) as pool2:
                signal = pool2.map(
                    partial(
                        extract_channel_data, edf_file=edf_file, chn_list=chn_labels
                    ),
                    chn_lists,
                )

        # Create headers:
        with Pool(processes=processes) as pool2:
            headers = pool2.map(
                partial(
                    extract_channel_header,
                    original_headers=headers_orig,
                    chn_list=chn_labels,
                    channels_labels=labels,
                ),
                chn_lists,
            )

        # Edit headers to make them compliant with edf files
        for header in headers:
            if new_srate != None:
                # f.datarecord_duration gives the value is sec and setDatarecordDuration receives it in units
                # of 10 ms. Therefore: setDatarecordDuration = datarecord_duration*10^6 / 10
                # This actually is used to set the sample frequency as the max number that can be written in
                # headers is 254 (uses int8)
                datarecord_dur, edf_srate = get_srate_params(new_srate)
                edf_out.setDatarecordDuration(int(datarecord_dur * 100000))
                header["sample_rate"] = edf_srate
                header["sample_frequency"] = edf_srate
            header["physical_max"] = int(header["physical_max"])
            header["physical_min"] = int(header["physical_min"])
            if len(str(header["physical_max"])) > 8:
                header["physical_max"] = int(str(header["physical_max"])[0:8])
            if len(str(header["physical_min"])) > 8:
                header["physical_min"] = int(str(header["physical_min"])[0:8])

        edf_out.setSignalHeaders(headers)
        edf_out.writeSamples(signal)
        edf_out.close()
    except Exception:
        traceback.print_exc()
        edf_out.close()
        edf_in.close()


# Function to look for timestamps
def extract_time_ids(epoch_id, time_vector, timestamps_array, srate):
    temp = np.asfortranarray(np.subtract(time_vector, timestamps_array[epoch_id]))
    t_init_id = np.abs(temp).argmin()
    t_end_id = int(np.floor(t_init_id + 240 * srate + 1))  # 4 min = 240 s
    return (t_init_id, t_end_id)


# Function to extract epochs
def extract_channel_epoch(chn_number, edf_file, srate_data, time_ids):
    edf_in = pyedflib.EdfReader(edf_file)
    n_val = time_ids[1] - time_ids[0]
    signal = edf_in.readSignal(chn_number, start=time_ids[0], n=n_val)
    # Deallocate space in memory
    edf_in.close()
    return signal


# Function to create epochs and EDF file from them
def create_epoch_EDF(edf_file, timestamp, out_path, processes):
    try:
        edf_in = pyedflib.EdfReader(edf_file)
        # First import labels
        labels = edf_in.getSignalLabels()
        # Create file:
        edf_out = pyedflib.EdfWriter(
            out_path, len(labels), file_type=pyedflib.FILETYPE_EDFPLUS
        )
        # First set the data from the header of the edf file:
        edf_out.setHeader(edf_in.getHeader())
        # f.datarecord_duration gives the value is sec and setDatarecordDuration receives it in units
        # of 10 ms. Therefore: setDatarecordDuration = datarecord_duration*10^6 / 10
        edf_out.setDatarecordDuration(
            int(edf_in.datarecord_duration * 100000)
        )  # This actually is used to set the sample frequency
        # Set each channel info:
        # Sampling rate:
        srate = edf_in.getSampleFrequencies()[0] / edf_in.datarecord_duration
        # Build epochs
        N = edf_in.getNSamples()[0]
        # Time vector:
        t = np.arange(0, N + 1) / srate
        # The plus 1 is required since the t[id] indicates start points.
        # Time ids
        t_init_id = np.abs(np.subtract(t, timestamp)).argmin()
        t_end_id = int(np.floor(t_init_id + 240 * srate))  # TODO: customizable
        n_val = t_end_id - t_end_id
        t_end_id = int(t_end_id - n_val % edf_in.getSampleFrequencies()[0])
        t_ids = (t_init_id, t_end_id)
        # Relative initial time for epoch
        t_0 = t[np.abs(np.subtract(t, timestamp)).argmin()]
        edf_out.writeAnnotation(0, -1, "Recording starts")
        # Headers for unipolar case
        headers = edf_in.getSignalHeaders()
        # Close file
        edf_in.close()
        # Extract channel information:
        chn_lists = range(len(labels))
        # Unipolar case:
        with Pool(processes=processes) as pool2:
            channel_data = pool2.map(
                partial(
                    extract_channel_epoch,
                    edf_file=edf_file,
                    srate_data=srate,
                    time_ids=t_ids,
                ),
                chn_lists,
            )
        # Edit headers to make them compliant with edf files
        for header in headers:
            header["physical_max"] = int(header["physical_max"])
            header["physical_min"] = int(header["physical_min"])
            if len(str(header["physical_max"])) > 8:
                header["physical_max"] = int(str(header["physical_max"])[0:8])
            if len(str(header["physical_min"])) > 8:
                header["physical_min"] = int(str(header["physical_min"])[0:8])

        edf_out.setSignalHeaders(headers)
        edf_out.writeSamples(channel_data)
        # Write annotations
        edf_out.writeAnnotation(t[t_init_id] - t_0, -1, f"Epoch starts.")
        edf_out.writeAnnotation(t[t_end_id + 1] - t_0, -1, f"Epoch ends.")
        # The plus 1 is required since the t[id] indicates start points, for the end epoch,
        # we want the final of the last datarecord!
        # Deallocate space in memory
        del t
        edf_out.close()
    except Exception:
        traceback.print_exc()
        edf_out.close()
        edf_in.close()


# # Function to extract epochs
# def extract_channel_epoch2(chn_number, edf_file, timestamp):
#     edf_in = EDFreader(edf_file)
#     # Sampling rate: based on channel
#     srate = edf_in.getSampleFrequency(chn_number)
#     # Build epochs
#     N = edf_in.getTotalSamples(chn_number)
#     # Time vector:
#     t = np.arange(0, N)/srate
#     # Time ids
#     t_init_id = np.abs(np.subtract(t,timestamp)).argmin()
#     t_end_id = int(np.floor(t_init_id+240*srate+1)) # TODO: customizable
#     n_val = t_end_id-t_init_id
#     # Read signal
#     dbuf = np.empty(n_val, dtype = np.float_)
#     edf_in.fseek(1, t_init_id, EDFreader.EDFSEEK_SET)
#     signal = edf_in.readSamples(1, dbuf, n_val)
#     # Deallocate space in memory
#     edf_in.close()
#     return signal

# # Function to create epochs and EDF file from them
# def create_epoch_EDF2(edf_file, timestamp, out_path, processes):
#     try:
#         edf_in = EDFreader(edf_file)
#         # First number of labels
#         n_labels = edf_in.getNumSignals()
#         # Create file:
#         edf_out = EDFwriter(out_path, EDFwriter.EDFLIB_FILETYPE_EDFPLUS, n_labels)
#         # Set headers
#         if edf_out.setStartDateTime(edf_in.getStartDateYear(), edf_in.getStartDateMonth(),
#                                     edf_in.getStartDateDay(), edf_in.getStartTimeHour(),
#                                     edf_in.getStartTimeMinute(), edf_in.getStartTimeSecond(),
#                                     edf_in.getStartTimeSubSecond()) != 0:
#             print("setStartDateTime() returned an error")
#             print(edf_in.getStartDateYear(), edf_in.getStartDateMonth(),
#                                     edf_in.getStartDateDay(), edf_in.getStartTimeHour(),
#                                     edf_in.getStartTimeMinute(), edf_in.getStartTimeSecond(),
#                                     edf_in.getStartTimeSubSecond()/1000)
#         if edf_out.setPatientCode(edf_in.getPatientCode()) != 0:
#             print("setPatientCode() returned an error")
#         edf_in_sex = edf_in.getPatientGender()
#         if edf_in_sex == 'Male':
#             if edf_out.setPatientGender(1) != 0:
#                 print("setPatientGender() returned an error")
#         elif edf_in_sex == 'Female':
#             if edf_out.setPatientGender(0) != 0:
#                 print("setPatientGender() returned an error")
#         else:
#             if edf_out.setPatientGender(2) != 0:
#                 print("setPatientGender() returned an error")
#         if edf_out.setPatientName(edf_in.getPatientName()) != 0:
#             print("setPatientName() returned an error")
#         if edf_out.setAdditionalPatientInfo(edf_in.getPatientAdditional()) != 0:
#             print("setAdditionalPatientInfo() returned an error")
#         if edf_out.setAdministrationCode(edf_in.getAdministrationCode()) != 0:
#             print("setAdministrationCode() returned an error")
#         if edf_out.setTechnician(edf_in.getTechnician()) != 0:
#             print("setTechnician() returned an error")
#         if edf_out.setEquipment(edf_in.getEquipment()) != 0:
#             print("setEquipment() returned an error")
#         if edf_out.setAdditionalRecordingInfo(edf_in.getRecordingAdditional()) != 0:
#             print("setAdditionalRecordingInfo() returned an error")
#         # Initial annotation
#         if edf_out.writeAnnotation(0, -1, "Recording starts") != 0:
#             print("writeAnnotation() returned an error")
#         # Close file to extract data
#         edf_in.close()
#         # Extract channel information:
#         chn_lists = range(n_labels)
#         with Pool(processes=processes) as pool2:
#             channel_data = pool2.map(partial(extract_channel_epoch2, edf_file=edf_file,
#                                             timestamp=timestamp), chn_lists)
#         # Write info channel based
#         # Reopen file
#         edf_in = EDFreader(edf_file)
#         for chan in chn_lists:
#             if edf_out.setPhysicalMaximum(chan, edf_in.getPhysicalMaximum(chan)) != 0:
#                 print("setPhysicalMaximum() returned an error")
#                 sys.exit()
#             if edf_out.setPhysicalMinimum(chan, edf_in.getPhysicalMinimum(chan)) != 0:
#                 print("setPhysicalMinimum() returned an error")
#                 sys.exit()
#             if edf_out.setDigitalMaximum(chan, edf_in.getDigitalMaximum(chan)) != 0:
#                 print("setDigitalMaximum() returned an error")
#                 sys.exit()
#             if edf_out.setDigitalMinimum(chan, edf_in.getDigitalMinimum(chan)) != 0:
#                 print("setDigitalMinimum() returned an error")
#                 sys.exit()
#             if edf_out.setPhysicalDimension(chan, edf_in.getPhysicalDimension(chan)) != 0:
#                 print("setPhysicalDimension() returned an error")
#                 sys.exit()
#             if edf_out.setSampleFrequency(chan, edf_in.getSampleFrequency(chan)) != 0:
#                 print("setSampleFrequency() returned an error")
#                 sys.exit()
#             if edf_out.setSignalLabel(chan, edf_in.getSignalLabel(chan)) != 0:
#                 print("setSignalLabel() returned an error")
#                 sys.exit()
#             if edf_out.setPreFilter(chan, edf_in.getPreFilter(chan)) != 0:
#                 print("setPreFilter() returned an error")
#                 sys.exit()
#             if edf_out.setTransducer(chan, edf_in.getTransducer(chan)) != 0:
#                 print("setTransducer() returned an error")
#                 sys.exit()
#             if edf_out.writeSamples(channel_data[chan]) != 0:
#                 print("writeSamples() returned error: %d" %(edf_out.writeSamples()))
#                 break

#         # Headers

#         # Deallocate space in memory
#         del t
#         edf_out.close()
#         edf_in.close()
#     except Exception:
#         traceback.print_exc()
#         edf_out.close()
#         edf_in.close()
