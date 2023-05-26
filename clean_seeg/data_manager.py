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
    # Only create bipolar channels with the closest or second closest electrode
    if abs(int(channels_dict[dict_key][id])-int(channels_dict[dict_key][id+1])) <= 2:
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
    # Only create bipolar channels with the closest or second closest electrode
    if abs(int(channels_dict[dict_key][id])-int(channels_dict[dict_key][id+1])) <= 2:
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
def create_bipolars(electrodes_df, electrodes_edf, processes, df_cols = None):
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
    channels = {}
    # Try to options of labels
    pattern1 = r'([A-Z0-9]+[-]+)(\d+)$' # for electrodes like 'LOpS-10' or 'LOpS1-10'
    pattern2 = r'([A-Z]+[-]*)(\d+)$' # for electrodes like 'LOpS10'
    # Extract channels info
    for electrode in electrodes_df[df_cols['label']].values:
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
    # Remove None elements
    bipolar_info_dicts = [element for element in bipolar_info_dicts if element is not None]
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
    if type(chn_list[0]) == list or type(chn_list[0]) == tuple: # Case 1: Bipolar
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
    elif type(chn_list[0]) == str: # Case 2: unipolar
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
    if type(chn_list[0]) == list or type(chn_list[0]) == tuple: # Case 1: Rereferencing to bipolar
        # Get indexes of channels
        chn1_id = channels_labels.index(chn_list[chn_number][1])
        # Update header
        chn_header = original_headers[chn1_id]
        chn_header['label'] = chn_list[chn_number][0]
    elif type(chn_list[0]) == str: # Case 2: no rereferencing
        # Get indexes of channel
        chn_id = channels_labels.index(chn_list[chn_number])
        # Copy header
        chn_header = original_headers[chn_id]
    return chn_header

# Function to get label map
def get_colors_labels():
    __location__ = os.path.realpath(
    os.path.join(os.getcwd(), os.path.dirname(__file__)))
    path_LUT = os.path.join(__location__, 'FreeSurferColorLUT.txt')
    with open(path_LUT, 'r') as f:
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
def get_chn_info(csv_file, electrodes_edf, df_cols = None): #, conf = 'unipolar'
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
    # print(list(df_cols.values()))
    important_data = df[list(df_cols.values())]
    important_data.reset_index() # make sure indexes pair with number of rows
    elec_df = pd.DataFrame(columns=list(df_cols.keys()))
    for index in range(len(important_data)):
        if important_data.loc[index, df_cols['label']] in electrodes_edf:
            tmp_df = pd.DataFrame([important_data.loc[index, list(df_cols.keys())].values], columns=list(df_cols.keys()))
            elec_df = pd.concat([elec_df, tmp_df], axis=0)
    del tmp_df
    # reset index
    elec_df = elec_df.reset_index(drop=True)
    # if conf == 'unipolar':
    chn_list = elec_df[df_cols['label']].values.tolist()
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
def get_srate_params(srate, max_srate=254): # max_srate=254 as uses int8
    # Convert sampling rate to integer
    srates = [srate, int(srate)] # try first on srate, if not, on the int
    for rate in srates:
        n_min = int(np.floor(rate/max_srate))
        for i in range(n_min+1, int(rate/2)):
            if rate % i == 0:
                duration = 1/i
                if len(str(duration).split(".")[1]) <= 4:
                    return 1/i, rate/i
    raise Exception('No appropriate parameters found in get_srate_params')
        
# Function to create EDF file based in channel list
def create_EDF(edf_file, out_path, processes, chn_labels=None, signal=None, n_removed = None, new_srate=None):
    import re
    # If signal == None, data is extracted from edf file based on the chn_labels
    try:
        edf_in = pyedflib.EdfReader(edf_file)
        # First import labels
        labels = edf_in.getSignalLabels()
        # Create file:
        if chn_labels==None:
            chn_labels = labels
        edf_out = pyedflib.EdfWriter(out_path, len(chn_labels), file_type=pyedflib.FILETYPE_EDFPLUS)
        # First set the data from the header of the edf file:
        edf_out.setHeader(edf_in.getHeader())
        headers_orig = edf_in.getSignalHeaders()

        # Set each channel info:
        if signal == None:
            N = edf_in.getNSamples()[0]
        else:
            N = len(signal[0]) # 0 index chosen randomly
        # Sampling rate:
        if new_srate == None:
            # f.datarecord_duration gives the value is sec and setDatarecordDuration receives it in units
            # of 10 ms. Therefore: setDatarecordDuration = datarecord_duration*10^6 / 10
            # This actually is used to set the sample frequency as the max number that can be written in headers is 254 (uses int8)
            edf_out.setDatarecordDuration(int(edf_in.datarecord_duration*100000)) 
            srate = edf_in.getSampleFrequencies()[0]/edf_in.datarecord_duration
        else:
            srate = new_srate
        
        # Annotations
        annot_orig = edf_in.readAnnotations()
        # Fill n_removed if not indicated
        # if n_removed == None:
        #     n_removed = np.zeros(len(annot_orig)+1)
        
        # Close file
        edf_in.close()
        print(f'len annot {len(annot_orig)}')
        # print(f'len n_rem {len(n_removed)}')
        # Write annotations
        t = np.arange(0, N)/srate
        # pattern_start = r'Epoch #\d starts.'
        # pattern_end = r'Epoch #\d ends.'
        # pattern_id = 0
        for annot_id in np.arange(len(annot_orig)):
            # print(annot_id)
            t_id = np.abs(np.subtract(t, annot_orig[0][annot_id])).argmin()
            t_annot = t[t_id]
            edf_out.writeAnnotation(t_annot, annot_orig[1][annot_id], annot_orig[2][annot_id])
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
        print('Channel part')
        if signal == None: # extracting data from edf file
            # Create bipolar signals:
            # print(chn_lists)
            with Pool(processes=processes) as pool2:
                signal = pool2.map(partial(extract_channel_data, edf_file=edf_file,
                                                chn_list=chn_labels), chn_lists)
        
        # Create headers:
        with Pool(processes=processes) as pool2:
            headers = pool2.map(partial(extract_channel_header, 
                                      original_headers=headers_orig,
                                      chn_list=chn_labels,
                                      channels_labels=labels), chn_lists)
        
        # Edit headers to make them compliant with edf files
        for header in headers:
            if new_srate != None:
                # f.datarecord_duration gives the value is sec and setDatarecordDuration receives it in units
                # of 10 ms. Therefore: setDatarecordDuration = datarecord_duration*10^6 / 10
                # This actually is used to set the sample frequency as the max number that can be written in 
                # headers is 254 (uses int8)
                datarecord_dur, edf_srate = get_srate_params(new_srate)
                edf_out.setDatarecordDuration(int(datarecord_dur*100000)) 
                header['sample_rate'] = edf_srate
                header['sample_frequency'] = edf_srate
            header['physical_max'] = int(header['physical_max'])
            header['physical_min'] = int(header['physical_min'])
            if len(str(header['physical_max']))>8:
                header['physical_max'] = int(str(header['physical_max'])[0:8])
            if len(str(header['physical_min']))>8:
                header['physical_min'] = int(str(header['physical_min'])[0:8])
        
        edf_out.setSignalHeaders(headers)
        edf_out.writeSamples(signal)
        edf_out.close()
    except Exception:
        traceback.print_exc()
        edf_out.close()
        edf_in.close()


# Function to look for timestamps
def extract_time_ids(epoch_id, time_vector, timestamps_array, srate):
    temp = np.asfortranarray(np.subtract(time_vector,timestamps_array[epoch_id]))
    t_init_id = np.abs(temp).argmin() 
    t_end_id = int(np.floor(t_init_id+240*srate+1)) # 4 min = 240 s
    return (t_init_id, t_end_id)

# Function to extract epochs
def extract_channel_epoch(chn_number, edf_file, srate_data, time_ids):
    edf_in = pyedflib.EdfReader(edf_file)
    n_val = time_ids[1]-time_ids[0]
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
        edf_out = pyedflib.EdfWriter(out_path, len(labels), file_type=pyedflib.FILETYPE_EDFPLUS)
        # First set the data from the header of the edf file:
        edf_out.setHeader(edf_in.getHeader())
        # f.datarecord_duration gives the value is sec and setDatarecordDuration receives it in units
        # of 10 ms. Therefore: setDatarecordDuration = datarecord_duration*10^6 / 10
        edf_out.setDatarecordDuration(int(edf_in.datarecord_duration*100000)) # This actually is used to set the sample frequency
        # Set each channel info:
        # Sampling rate:
        srate = edf_in.getSampleFrequencies()[0]/edf_in.datarecord_duration
        # Build epochs
        N = edf_in.getNSamples()[0]
        # Time vector:
        t = np.arange(0, N)/srate
        # Time ids
        t_init_id = np.abs(np.subtract(t,timestamp)).argmin()
        t_end_id = int(np.floor(t_init_id+240*srate+1)) # TODO: customizable
        t_ids = (t_init_id, t_end_id)
        # Relative initial time for epoch
        t_0 = t[np.abs(np.subtract(t,timestamp)).argmin()]
        edf_out.writeAnnotation(0, -1, "Recording starts")
        # Headers for unipolar case
        headers = edf_in.getSignalHeaders()
        # Close file
        edf_in.close()
        # Extract channel information:
        chn_lists = range(len(labels))
        # Unipolar case:
        with Pool(processes=processes) as pool2:
            channel_data = pool2.map(partial(extract_channel_epoch, edf_file=edf_file, 
                                            srate_data=srate, time_ids=t_ids), chn_lists)
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
        # Write annotations
        edf_out.writeAnnotation(t[t_init_id]-t_0, -1, f"Epoch starts.")
        edf_out.writeAnnotation(t[t_end_id]-t_0, -1, f"Epoch ends.")
        # Deallocate space in memory
        del t
        edf_out.close()
    except Exception:
        traceback.print_exc()
        edf_out.close()
        edf_in.close()