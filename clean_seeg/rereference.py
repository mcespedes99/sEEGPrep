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

# Function to create bipolar channels from given unipolars
def create_bipolar_comb(id, dict_key, channels_dict):
    bipolar_chn = dict_key+channels_dict[dict_key][id]+'-'+channels_dict[dict_key][id+1] 
    chn1_label = dict_key+channels_dict[dict_key][id]
    chn2_label = dict_key+channels_dict[dict_key][id+1]
    return (bipolar_chn, chn1_label, chn2_label)

def create_bipolar_combi(id, bipolar_list):
    return bipolar_list[id]

# Function to discard 
def apply_bipolar_criteria(df_bipolar, bipolar_list, processes):
    non_white_matter_unknown_bool = df_bipolar['Label'].str.contains('White-Matter|Unknown',
                                                                    case=False, regex=True)==False
    ids = df_bipolar.loc[non_white_matter_unknown_bool, 'Label'].index.values.tolist()
    # print(ids)
    # Extract bipolar list (original channels + bipolar channel)
    with Pool(processes=processes) as pool:
            filtered_list = pool.map(partial(create_bipolar_combi, bipolar_list=bipolar_list), 
                                    ids)
    return filtered_list

# Function to extract position of bipolar channels
def bipolar_info(id, dict_key, channels_dict, elec_pos, df_cols):
    # Bipolar channel label
    bipolar_chn = dict_key+channels_dict[dict_key][id]+'-'+channels_dict[dict_key][id+1] 
    # Channels' labels
    chn1_label = dict_key+channels_dict[dict_key][id]
    chn2_label = dict_key+channels_dict[dict_key][id+1]
    # Extract positions
    inf_chn1 = elec_pos.loc[elec_pos[df_cols['label']] == chn1_label]
    inf_chn2 = elec_pos.loc[elec_pos[df_cols['label']] == chn2_label]
    # print(inf_chn1)
    data = {
        'type': inf_chn1[df_cols['type']].values[0],
        'group': inf_chn1[df_cols['group']].values[0],
        'label': bipolar_chn,
        'x': (inf_chn1[df_cols['x']].values[0] + inf_chn2[df_cols['x']].values[0])/2,
        'y': (inf_chn1[df_cols['y']].values[0] + inf_chn2[df_cols['y']].values[0])/2,
        'z': (inf_chn1[df_cols['z']].values[0] + inf_chn2[df_cols['z']].values[0])/2
    }
    return data
    

# Function to create a bipolar channel list from 
def create_bipolars(electrodes_df, processes, df_cols = None):
    # df_cols (dict) = {type_record, label, x, y, z, group}
    # the dict can be in any order. The df will have the structure given
    # by the dict
    if df_cols == None:
        df_cols = {
            'type': 'type',
            'label': 'label',
            'x': 'x',
            'y': 'y',
            'z': 'z',
            'group': 'group',
        }
    # print(df_cols['group'])
    channels = dict((label,[]) for label in electrodes_df[df_cols['group']].unique())
    pattern = r'([A-Z]+)(\d+)'
    electrode_labels = electrodes_df[df_cols['label']].values
    # Extract channels info
    for electrode in electrodes_df[df_cols['label']].values:
        match = re.match(pattern, electrode, re.IGNORECASE)
        channels[match.group(1)].append(match.group(2))
    # print(channels)
    # Create new list
    bipolar_list = []
    bipolar_info_dicts = []
    for key in  channels.keys():
        # Extract bipolar list (original channels + bipolar channel)
        with Pool(processes=processes) as pool:
                bipolar_list = bipolar_list + pool.map(partial(create_bipolar_comb, dict_key=key, channels_dict=channels), 
                                        list(range(len(channels[key])-1)))
        # Extract bipolar channels information
        with Pool(processes=processes) as pool:
                bipolar_info_dicts = bipolar_info_dicts + pool.map(partial(bipolar_info, dict_key=key, channels_dict=channels,
                                                                           elec_pos=electrodes_df, df_cols=df_cols), 
                                                                   list(range(len(channels[key])-1)))
        bipolar_elec = pd.DataFrame(columns=  list(df_cols.keys()))
        data = pd.DataFrame(bipolar_info_dicts)
        bipolar_elec = pd.concat([bipolar_elec, data], ignore_index=True)
    return bipolar_list, bipolar_elec

# Function to extract info from each channel - bipolar
def extract_channel_data(chn_number, edf_file, srate_data, bipolar_list):
    edf_in = pyedflib.EdfReader(edf_file)
    # Get labels from original edf file
    channels_labels = edf_in.getSignalLabels()
    # Get indexes of channels
    # print('new iter')
    chn1_id = channels_labels.index(bipolar_list[chn_number][1])
    chn2_id = channels_labels.index(bipolar_list[chn_number][2])
    # print(chn_number)
    signal_chn1 = edf_in.readSignal(chn1_id)
    signal_chn2 = edf_in.readSignal(chn2_id)
    chn_data = signal_chn1 - signal_chn2
    # Deallocate space in memory
    edf_in.close()
    del signal_chn1
    del signal_chn2
    return chn_data


# Function to extract headers for bipolar channels
def extract_channel_header(chn_number, original_headers, bipolar_list, channels_labels):
    # Get indexes of channels
    chn1_id = channels_labels.index(bipolar_list[chn_number][1])
    # Update header
    chn_header = original_headers[chn1_id]
    chn_header['label'] = bipolar_list[chn_number][0]
    return chn_header

# Function to get label map
def get_colors_labels():
    with open('/home/mcesped/scratch/code/sEEGPrep/clean_seeg/FreeSurferColorLUT.txt', 'r') as f:
        raw_lut = f.readlines()

    # read and process line by line
    label_map = pd.DataFrame(columns=['Label', 'R', 'G', 'B'])
    for line in raw_lut:
        # Remove empty spaces
        line = line.strip()
        if not (line.startswith('#') or not line):
            s = line.split()
            # info = list(filter(None, info))
            id = int(s[0])
            info_s = {
                'Label': s[1],
                'R': int(s[2]),
                'G': int(s[3]),
                'B': int(s[4])
            }
            # info_s['A'] = 0 if (info_s['R']==0 & info_s['G']==0 & info_s['B']==0) else 255
            info_s = pd.DataFrame(info_s, index=[id])
            label_map = pd.concat([label_map,info_s], axis=0)
        label_map[['R','G','B']] = label_map[['R','G','B']].astype('int64')
    return label_map

# Function to get label id based on parcellation obj
def get_electrodes_id(parc, elec_df, non_cont_to_cont_tf, df_cols):
    # Load data of parcellations
    data_parc = np.asarray(parc.dataobj)
    # Coordinates in MRI RAS
    mri_ras_mm = elec_df[[df_cols['x'],df_cols['y'],df_cols['z']]].values
    # print(mri_ras_mm)
    # Transform from contrast mri ras to non-contrast MRI ras
    mri_ras_mm = mne.transforms.apply_trans(non_cont_to_cont_tf, mri_ras_mm)
    # print(mri_ras_mm)
    # To voxels
    inv_affine = np.linalg.inv(parc.affine)
    # here's where the interpolation should be performed!!
    vox = (mne.transforms.apply_trans(inv_affine, mri_ras_mm)).astype(int)
    # print(vox)
    id = data_parc[vox[:,0], vox[:,1], vox[:,2]]
    return id

# Function to get rgb values for each contact
def get_label_rgb(parc, elec_df, non_cont_to_cont_tf, label_map, df_cols):
    # vox, data_parc = ras2vox(parc, elec_df, non_cont_to_cont_tf)
    id = get_electrodes_id(parc, elec_df, non_cont_to_cont_tf, df_cols)
    vals = label_map.loc[id, ['Label','R', 'G', 'B']].to_numpy()
    vals = np.c_[id,vals]
    vals = pd.DataFrame(data=vals, columns=['Label ID','Label','R', 'G', 'B'])
    vals = pd.concat([elec_df, vals], axis=1)
    return vals

# Function to read matrix
def readRegMatrix(trsfPath):
	with open(trsfPath) as (f):
		return np.loadtxt(f.readlines())

# Function to create tsv with bipolar channels info
def extract_location(parc_path, noncon_to_con_tf_path, chn_info_df, df_cols):
    import os
    # Load labels from LUT file
    labels = get_colors_labels()
    # Load parcellation file
    parc_obj = nb.load(parc_path)
    data_parc = np.asarray(parc_obj.dataobj)
    # The transform file goes from contrast to non-contrast. The tfm, when loaded in slicer actually goes
    # from non-contrast to contrast but the txt is inversed!
    t1_transform=readRegMatrix(noncon_to_con_tf_path)
    # Create df
    df = get_label_rgb(parc_obj, chn_info_df, t1_transform, labels, df_cols)
    # if not os.path.exists(out_tsv_name):
    #     df.to_csv(out_tsv_name, sep = '\t')
    return df

# Function to extract useful information from csv file
def get_chn_info(csv_file, df_cols = None):
    df = pd.read_csv(csv_file, sep='\t')
    # df_cols (dict) = {type_record, label, x, y, z, group}
    # the dict can be in any order. The df will have the structure given
    # by the dict
    if df_cols == None:
        df_cols = {
            'type': 'type',
            'label': 'label',
            'x': 'x',
            'y': 'y',
            'z': 'z',
            'group': 'group',
        }
    print(list(df_cols.values()))
    important_data = df[list(df_cols.values())].values
    elec_df = pd.DataFrame(columns=  list(df_cols.keys()), data = important_data)
    chn_list = df[df_cols['label']].values.tolist()
    return elec_df, chn_list, df_cols
    

# Function to create EDF file with bipolar data
def create_EDF(edf_file, bipolar_channels, out_path, processes):
    try:
        edf_in = pyedflib.EdfReader(edf_file)
        # First import labels
        labels = edf_in.getSignalLabels()
        # Create file:
        edf_out = pyedflib.EdfWriter(out_path, len(bipolar_channels), file_type=pyedflib.FILETYPE_EDFPLUS)
        # First set the data from the header of the edf file:
        edf_out.setHeader(edf_in.getHeader())
        headers_orig = edf_in.getSignalHeaders()
        # f.datarecord_duration gives the value is sec and setDatarecordDuration receives it in units
        # of 10 ms. Therefore: setDatarecordDuration = datarecord_duration*10^6 / 10
        edf_out.setDatarecordDuration(int(edf_in.datarecord_duration*100000)) # This actually is used to set the sample frequency
        # Set each channel info:
        # Sampling rate:
        srate = edf_in.getSampleFrequencies()[0]/edf_in.datarecord_duration
        # Build epochs
        N = edf_in.getNSamples()[0]
        
        # Annotations
        annot_orig = edf_in.readAnnotations()
        
        # Close file
        edf_in.close()
        # Write annotations
        for annot_id in np.arange(len(annot_orig)):
            # print(annot_id)
            edf_out.writeAnnotation(max(0, annot_orig[0][annot_id]), annot_orig[1][annot_id], annot_orig[2][annot_id])
        
        # Extract channel information:
        chn_lists = range(len(bipolar_channels)) #for bipolar
        print('Channel part')
        # Create bipolar signals:
        print(chn_lists)
        with Pool(processes=processes) as pool2:
            channel_data = pool2.map(partial(extract_channel_data, edf_file=edf_file, 
                                            srate_data=srate,
                                            bipolar_list=bipolar_channels), chn_lists)
        # Create headers:
        with Pool(processes=processes) as pool2:
            headers = pool2.map(partial(extract_channel_header, 
                                      original_headers=headers_orig,
                                      bipolar_list=bipolar_channels,
                                      channels_labels=labels), chn_lists)
        
        # Edit headers to make them compliant with edf files
        for header in headers:
            header['physical_max'] = int(header['physical_max'])
            header['physical_min'] = int(header['physical_min'])
            if len(str(header['physical_max']))>8:
                header['physical_max'] = int(str(header['physical_max'])[0:8])
            if len(str(header['physical_min']))>8:
                header['physical_min'] = int(str(header['physical_min'])[0:8])
        
        edf_out.setSignalHeaders(headers)
        edf_out.writeSamples(channel_data)
        edf_out.close()
    except Exception:
        traceback.print_exc()
        edf_out.close()
        edf_in.close()
