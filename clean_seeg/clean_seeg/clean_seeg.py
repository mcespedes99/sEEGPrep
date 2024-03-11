from .utils import (
    get_montage,
    get_chn_labels,
    plot_filter,
    get_chn_positions,
    downsampling,
    remove_trend,
)
from .clean_autoreject import create_mne_epochs, run_autoreject
from .clean_flatlines import clean_flatlines
from .clean_PLI import removePLI_chns, zapline, cleanline, notch_filt
from .data_manager import (
    create_EDF,
    create_bipolars,
    extract_location,
    apply_bipolar_criteria,
    get_chn_info,
    create_epoch_EDF,
)
import logging
import pandas as pd


class cleanSEEG:
    def __init__(
        self,
        edf_path,
        RmTrendMethod="HighPass",
        methodPLI="Zapline",
        lineFreq=60,
        bandwidth=4,
        n_harmonics=3,
        noiseDetect=True,
        highpass=[0.5, 1],
        tfm=[],
        epoch_autoreject=5,  # Epoch length for autoreject
        processes=None,
    ):
        import pyedflib

        self.edf_path = edf_path
        self.RmTrendMethod = RmTrendMethod
        self.methodPLI = methodPLI
        self.lineFreq = lineFreq
        self.bandwidth = bandwidth
        self.n_harmonics = n_harmonics
        self.noiseDetect = noiseDetect
        self.highpass = highpass  # Set to None to shut down
        self.tfm = tfm  # has to be a list of tuple with format: (path, invert), with 'path' as str or tfm and invert as boolean
        self.epoch_autoreject = epoch_autoreject
        self.processes = processes
        # Extra params
        self.reref_chn_list = []
        self.reref_chn_df = []
        self.inter_edf = None
        self.clean_edf = None
        self.subject = (None,)
        self.subjects_dir = None

        # Find sample rate
        edf_in = pyedflib.EdfReader(edf_path)
        self.srate = edf_in.getSampleFrequencies()[0] / edf_in.datarecord_duration
        edf_in.close()

    # Epoch extraction
    def extract_epochs(
        self,
        n_samples=None,
        event_dur=None,
        event_start=None,
        out_root=None,
        out_files=None,
        tmpdir=None,
        snakemake=False,
        return_report=True
    ):
        import pyedflib
        import shutil
        import os
        import re
        import snakebids
        import bids
        import numpy as np

        assert n_samples or event_dur
        # json obj with report info
        report = {
            'Epoch start event': [],
            'Relative start time': [],
            'Duration': [],
            'Clip number' : []
        }
        # Find indexes from events
        f = pyedflib.EdfReader(self.edf_path)
        # Number of samples in edf
        N = f.getNSamples()[0]
        # Get number of samples if not defined
        # Sampling rate:
        srate = f.getSampleFrequencies()[0] / f.datarecord_duration
        if not n_samples:
            # Get number of samples based on event duration
            n_samples = int(event_dur * srate)
        # Check than n_samples is not bigger than N
        if n_samples > N:
            print('Number of samples bigger than number of points in file, changing to the max.\n', flush=True)
            n_samples = N
        # Adjust n_samples depending on the datarecord duration and the number of samples in file
        remainder = n_samples % f.getSampleFrequencies()[0]
        if remainder != 0:
            n_samples += f.getSampleFrequencies()[0] - remainder
        
        # Find timestamps
        t = np.arange(0, N) / srate
        if event_start is not None:
            # Two possibilities: label or index
            if isinstance(event_start, int):
                time_stamps_init = t[event_start]
            else:    
                id = [
                    value[0]
                    for value in enumerate(f.readAnnotations()[2])
                    if re.match(event_start, value[1], re.IGNORECASE)
                ]
                # Create df with annotations
                onset_list = f.readAnnotations()[0]
                # Find time stamps where the 'awake trigger' event is happening
                time_stamps_init = onset_list[id]
            # Make sure it's a list
            if not hasattr(time_stamps_init, '__iter__'):
                time_stamps_init = [time_stamps_init]
        else:
            time_stamps_init = t[::n_samples]

        f.close()
        # print(time_stamps)
        # Copy file to local scratch if possible
        new_edf = None
        if tmpdir != None and os.path.exists(tmpdir):
            # Copy file to local scratch
            print("here")
            file_name = os.path.basename(self.edf_path)
            new_edf = os.path.join(tmpdir, file_name)
            shutil.copy(self.edf_path, new_edf)

        # Extract entities from input path
        entities = bids.layout.parse_file_entities(self.edf_path)
        # Combine 'extension' with 'suffix' and delete the first one
        if snakemake:
            entities["suffix"] = entities["suffix"] +'_tmp'+ entities["extension"]
        else:
            entities["suffix"] = entities["suffix"] + entities["extension"]
        del entities["extension"]
        # Add 'task'
        entities["rec"] = "clip"
        # Here call function to create new EDF file
        for index, event_timestamp in enumerate(time_stamps_init):
            if out_files is None:
                # New file name
                entities["clip"] = f"{index+1:02}"
                out_edf_path = os.path.basename(
                    snakebids.bids(root=out_root, **entities)
                )
                out_edf_path = os.path.join(out_root, out_edf_path)
                print(out_edf_path)
            else:
                # Filter based on regex
                reg = re.compile(f"clip-{index+1:02}")
                out_edf_path = list(filter(reg.search, out_files))[0]
            if new_edf == None:
                time_init, duration = create_epoch_EDF(
                                        self.edf_path,
                                        event_timestamp,
                                        n_samples,
                                        out_edf_path,
                                        self.processes,
                )
            else:
                print("aqui")
                time_init, duration = create_epoch_EDF(
                    new_edf, event_timestamp, n_samples, out_edf_path, self.processes
                )
            # Save report
            report["Epoch start event"].append(event_start)
            report['Relative start time'].append(time_init)
            report["Duration"].append(duration)
            report['Clip number'].append(f"clip-{index+1:02}")
        if new_edf != None:
            print("delete")
            os.remove(new_edf)
        if return_report:
            return pd.DataFrame(report)

    # Reference function
    def rereference(
        self,
        electrodes_tsv,
        out_edf_path,
        write_tsv=False,
        out_tsv_path=None,
        return_report=True,
    ):
        # TODO: include use_clean feature!
        import os
        import pyedflib
        import numpy as np

        # Manage some errors:
        if write_tsv and out_tsv_path == None:
            raise Exception("To write the tsv file, please indicate an output path.")
        if electrodes_tsv == None:
            raise Exception(
                "Please indicate the path to a tsv file with channels information."
            )
        # TODO: separate rerefering from label map creation
        print("Running rereference")
        # Extract info about electrodes positions
        elec_pos = pd.read_csv(electrodes_tsv, sep="\t")
        # Extract channels present in edf file (might not be all)
        edf = pyedflib.EdfReader(self.edf_path)
        elec_edf = edf.getSignalLabels()
        edf.close()
        # First check with electrodes are in electrodes.tsv and EDF
        _, discarded_labels = get_chn_labels(electrodes_tsv, elec_edf)
        assert len(discarded_labels)<len(elec_pos) # If not, it means the edf if wrong format or it's already bipolar
        # Create bipolar combinations
        bipolar_channels, bipolar_info_df = create_bipolars(
            elec_pos, elec_edf, self.processes
        )
        # Save values in the class
        self.reref_chn_list = bipolar_channels
        self.reref_chn_df = bipolar_info_df
        # Write tsv
        if write_tsv:  # and not os.path.exists(out_tsv_name):
            if os.path.exists(out_tsv_path):
                logging.warning(f"tsv file {out_tsv_path} will be overwritten.")
            bipolar_info_df.to_csv(out_tsv_path, index=False, sep="\t")

        # print(bipolar_channels)

        # Create new EDF file
        if os.path.exists(out_edf_path):
            logging.warning(f"edf file {out_edf_path} will be overwritten.")
        create_EDF(
            self.edf_path, out_edf_path, self.processes, chn_labels=bipolar_channels
        )
        self.rereference_edf = out_edf_path
        if return_report:
            bipolar_combs = {
                'Bipolar channel': [values[0] for values in bipolar_channels],
                'Unipolar channels': [values[1:] for values in bipolar_channels]
            }
            report = {
                'Discarded channels': discarded_labels
            }
            return pd.DataFrame(bipolar_combs), report

    def identify_regions(
        self,
        electrodes_tsv,
        aparc_aseg_path,
        colortable_file,  # colortable has to be a tsv with at least index, name
        reference, # 
        use_reref=True,
        write_tsv=False,
        out_tsv_path=None,
        discard_wm_un=False,
        write_edf=False,
        out_edf_path=None,
        vol_version=False,
        json_out=None,
        mask_out=None, # to save mask
        return_report=True
    ):
        import os
        import pyedflib
        import json

        assert reference in ['bipolar', 'unipolar']

        # Manage exception
        if electrodes_tsv == None:
            raise Exception(
                "Please indicate the path to a tsv file with channels information."
            )
        # discard_wm_un discards white matter and unknown from edf file (not from csv)
        # Manage a few exceptions
        if write_edf and out_edf_path == None:
            raise Exception(
                "Please indicate a value for out_edf_path or set write_edf to False."
            )

        if use_reref:
            if len(self.reref_chn_list) == 0 or len(self.reref_chn_df) == 0:
                raise Exception(
                    "Please run rereference first or set use_reref to False."
                )
            chn_info_df = self.reref_chn_df
            chn_list = self.reref_chn_list
        # Extract electrodes information if not using rereference results
        else:
            logging.info(f"Extracting channel positions from {electrodes_tsv}")
            # Extract info about electrodes positions
            elec_pos = pd.read_csv(electrodes_tsv, sep="\t")
            # Filter any channels with nan vals
            nan_rows = elec_pos[['x', 'y', 'z']].isna().any(axis=1)
            elec_pos = elec_pos.loc[~nan_rows]
            dropped_chns = elec_pos.loc[nan_rows, 'name'].values.tolist()
            # First extract list of channels present in edf file
            edf = pyedflib.EdfReader(self.edf_path)
            elec_edf = edf.getSignalLabels()
            edf.close()
            # Check if needs to preprocess electrodes.tsv file
            if reference=='bipolar':
                # Needs to convert electrodes df to bipolar
                # Create bipolar combinations
                _, elec_pos = create_bipolars(
                    elec_pos, elec_edf, self.processes, compare_edf=False
                )

            chn_info_df, chn_list, discarded_chns = get_chn_info(
                elec_pos, elec_edf, reference
            )  # , conf=conf
        # Create tsv file with information about location of the channels
        df_location, regions_per_chn = extract_location(
            aparc_aseg_path,
            chn_info_df,
            self.tfm,
            colortable_file,
            reference,
            vol_version,
            mask_out
        )
        if write_tsv:
            if os.path.exists(out_tsv_path):
                logging.warning(f" tsv file {out_tsv_path} will be overwritten.")
            df_location.to_csv(out_tsv_path, index=False, sep="\t")
        # Write json
        if json_out and vol_version:
            with open(json_out, "w") as outfile:
                json.dump(regions_per_chn, outfile)
        # Discard data from white matter
        if discard_wm_un:
            chn_list = apply_bipolar_criteria(df_location, chn_list, self.processes)
            # print(chn_list)
            # Overwrite reref info if previously run
            if use_reref:
                self.reref_chn_list = chn_list
            # else:
            #     logging.critical('Tool currently does not write the edf file for unipolar cases.')
            #     return df_location

            # Write edf file if requested (only if the chn list changed!)
            # NOT WORKING FOR UNIPOLAR CASES! create_EDF only works for bipolar
            if write_edf:
                if os.path.exists(out_edf_path):
                    logging.warning(f"EDF file {out_edf_path} will be overwritten.")
                create_EDF(
                    self.edf_path, out_edf_path, self.processes, chn_labels=chn_list
                )
                self.rereference_edf = out_edf_path
        elif (not discard_wm_un) and write_edf:
            logging.warning(
                f"EDF file {out_edf_path} will only be copied as no updates have been made."
            )
        
        if return_report:
            unique_regions = df_location['region name'].unique()
            region_maps = {'Region': unique_regions, 'Number of channels': [], 'Channels':[]}
            for region in unique_regions:
                chns = df_location.loc[df_location['region name']==region, "name"].values
                region_maps['Number of channels'].append(len(chns))
                region_maps['Channels'].append(chns)
            df_regions = pd.DataFrame(region_maps)
            report = {
                'Discarded channels': discarded_chns+dropped_chns
            }
            return df_location, df_regions, report
        return df_location
        

    def downsample(
        self, channels_tsv, target_srate, write_edf=False, out_edf_path=None, return_report=True
    ):
        import pyedflib
        from multiprocessing.pool import Pool
        from multiprocessing import get_context
        from functools import partial
        import numpy as np
        import os

        # print('aqui')
        # print(self.edf_path)
        # Open edf file
        edf_in = pyedflib.EdfReader(self.edf_path)
        # Extract labels to only process signals in tsv file
        labels = edf_in.getSignalLabels()
        chn_labels, discarded_labels = get_chn_labels(channels_tsv, labels)
        # Number of samples
        N = edf_in.getNSamples()[0]
        # Create signal list
        channels = []
        # Extract signal per channel
        for chan in chn_labels:
            id_ch = labels.index(chan)
            channels.append(id_ch)
        edf_in.close()
        # create a process context. Refer to:
        # https://github.com/dask/dask/issues/3759
        ctx = get_context("spawn")
        with Pool(processes=self.processes, context=ctx) as pool:
            data_list, newSrate = zip(
                *pool.map(
                    partial(
                        downsampling,
                        edf_file=self.edf_path,
                        orig_srate=self.srate,
                        target_srate=target_srate,
                    ),
                    channels,
                )
            )
        # Reformatting to array
        n_samples = len(data_list[0])
        data_dnsampled = np.zeros((len(chn_labels), n_samples))
        for ch in np.arange(len(data_list)):
            data_dnsampled[ch, :] = data_list[ch]

        # Write edf if requested
        if write_edf:
            if os.path.exists(out_edf_path):
                logging.warning(f"EDF file {out_edf_path} will be overwritten.")
            create_EDF(
                self.edf_path,
                out_edf_path,
                self.processes,
                signal=data_list,
                new_srate=newSrate[0],
                chn_labels=chn_labels
            )
            del data_list
            self.downsample_edf = out_edf_path
        report = {
            'Original sampling rate': self.srate,
            'Target sampling rate': target_srate,
            'New sampling rate': newSrate[0],
            'Discarded channels':discarded_labels
        }
        if return_report:
            return data_dnsampled, pd.DataFrame([report])
        else:
            return data_dnsampled

    # Function to reject PLI
    def reject_PLI(self, channels_tsv, write_edf=False, out_edf_path=None, return_report=True):
        import pyedflib
        import numpy as np
        import traceback
        import os

        # Manage a few exceptions:
        if write_edf and out_edf_path == None:
            raise Exception(
                "EDF file with clean signal cannot be written without and appropiate out path"
            )
        # Open edf file
        edf_in = pyedflib.EdfReader(self.edf_path)
        # Extract labels
        labels = edf_in.getSignalLabels()
        chn_labels, discarded_labels = get_chn_labels(channels_tsv, labels)
        try:
            # Number of samples
            N = edf_in.getNSamples()[0]
            # Create signal list
            signal = []
            # Extract signal per channel
            for chan in chn_labels:
                id_ch = labels.index(chan)
                chn_sig = edf_in.readSignal(id_ch)
                signal.append(chn_sig)
            edf_in.close()
            # Convert signal to array
            signal = np.vstack(signal)
            # Remove line noise
            print("Removing line noise")
            clean = self.wrapper_PLI(signal)
            print("PLI removal completed.")
            # Write edfs if requested
            if write_edf:
                # Clean signal
                # convert to list
                clean_list = [clean[i, :] for i in range(clean.shape[0])]
                if os.path.exists(out_edf_path):
                    logging.warning(f"EDF file {out_edf_path} will be overwritten.")
                create_EDF(
                    self.edf_path,
                    out_edf_path,
                    self.processes,
                    chn_labels=chn_labels,
                    signal=clean_list,
                )
                del clean_list
                self.clean_edf = out_edf_path
            if return_report:
                report = dict()
                report['Method'] = self.methodPLI
                report['Line frequency (Hz)'] = self.lineFreq
                if self.methodPLI == "Cleanline":
                    report['Bandwidth (Hz)'] = self.bandwidth
                elif self.methodPLI == "NotchFilter":  # TODO: add bandwidth param
                    report['Number of Harmonics'] = self.n_harmonics
                elif self.methodPLI == "PLIremoval":
                    report['Number of Harmonics'] = 3
                    report['B0, Initial notch bandwidth of the frequency estimator'] = 100
                    report['Binf, Asymptotic notch bandwidth of the frequency estimator'] = 0.01
                    report['Bst, Rate of convergence to 95# of the asymptotic bandwidth Binf'] = 4
                    report['P0, Initial settling time of the frequency estimator'] = 0.1
                    report['Asymptotic settling time of the frequency estimator'] = 2
                    report['Pst, Rate of convergence to 95# of the asymptotic settling time'] = 5
                    report['W, Settling time of the amplitude and phase estimator'] = 2
                    # Hardcoded for now
                return clean, report
            return clean
        except:
            edf_in.close()
            print(traceback.format_exc())
            raise Exception

    # Wrapper to manage different PLI methods
    def wrapper_PLI(self, signal):
        if self.methodPLI == "Cleanline":
            signal = cleanline(
                signal.T, self.srate, processes=self.processes, bandwidth=self.bandwidth
            )
        elif self.methodPLI == "Zapline":
            signal = zapline(signal.T, self.lineFreq / self.srate, self.srate)
        elif self.methodPLI == "NotchFilter":  # TODO: add bandwidth param
            signal = notch_filt(
                signal.T, self.lineFreq, self.srate, n_harmonics=self.n_harmonics
            )
        elif self.methodPLI == "PLIremoval":
            signal = removePLI_chns(
                signal.T,
                self.srate,
                3,
                [100, 0.01, 4],
                [0.1, 2, 5],
                2,
                processes=self.processes,
                f_ac=self.lineFreq,
            )  # Hardcoded for now
        else:
            raise Exception("PLI method not valid.")
        signal = signal.T
        return signal

    def clean_epochs(
        self,
        channels_tsv,
        subject=None,
        subjects_dir=None,
        return_interpolated=False,
        write_edf_clean=False,
        out_edf_path_clean=None,
        write_tsv=False,
        out_tsv_path=None,
        write_edf_int=False,  # Not working for now
        out_edf_path_int=None,
        verbose=False,
    ):
        import pyedflib
        import numpy as np
        import pandas as pd
        import traceback
        import logging
        import os

        # Include attributes:
        self.subject = subject
        self.subjects_dir = subjects_dir
        # Manage a few exceptions:
        if write_edf_clean and out_edf_path_clean == None:
            raise Exception(
                "EDF file with clean signal cannot be written without and appropiate out path"
            )
        if write_edf_int and out_edf_path_int == None:
            raise Exception(
                "EDF file with interpolated signal cannot be written without and appropiate out path"
            )
        if write_tsv and out_tsv_path == None:
            raise Exception(
                "TSV file with noise information cannot be written without and appropiate out path"
            )
        if channels_tsv == None:
            raise Exception(
                "Please indicate the path to a tsv file with channels information."
            )  # Open edf file
        edf_in = pyedflib.EdfReader(self.edf_path)
        # Extract labels
        elec_edf = edf_in.getSignalLabels()
        # Begin by getting the position of the electrodes in RAS space
        chn_labels, discarded_labels = get_chn_labels(channels_tsv, elec_edf)
        try:
            # Number of samples
            N = edf_in.getNSamples()[0]
            ### THIS MUST BE CHANGED TO A MORE STANDARD APPROACH
            # Define length of epochs based on the first one
            # t_init_id = np.argmin((np.abs(t-t_init)))
            # t_end_id = np.argmin((np.abs(t-t_end)))
            # samples_epoch = t_end_id-t_init_id
            # print(samples_epoch)
            # Initiate csv epoch file
            cols = ["Epoch #", "Start ID", "End ID"] + chn_labels
            df_epochs = pd.DataFrame(columns=cols)
            # Count of removed elements per epoch
            n_removed = [0]
            # Initialize in 0 as per epoch, the start id will be n_orig[id]-n_removed[id] and the end id: n_orig[id]-n_removed[id+1]
            # Example: if 10 elements were removed in the first epoch and 12 in the 2nd, then the start id of the second should be
            # the original one minus 10 but the end should be the original minus 22
            # Create signal list
            signal = []
            # Extract signal per channel
            for chan in chn_labels:
                id_ch = elec_edf.index(chan)
                chn_sig = edf_in.readSignal(id_ch)
                signal.append(chn_sig)
            edf_in.close()
            # Convert signal to array
            signal = np.vstack(signal)
            # Run cleaning
            # First remove line noise
            # print('Removing line noise')
            # signal = self.wrapper_PLI(signal)
            # print('PLI removal completed.')

            # Second, identify noisy segments and highpass filter the data
            if self.noiseDetect and return_interpolated:
                clean, interpolated, df_epochs, filt_params = self.noiseDetect_raw(
                    signal, elec_edf, channels_tsv, return_interpolated=True, verbose=verbose
                )
                # Update number of removed elements after autoreject
                n_removed.append(clean.shape[-1] - interpolated.shape[-1])
            elif self.noiseDetect:
                clean, df_epochs, filt_params = self.noiseDetect_raw(
                    signal, elec_edf, channels_tsv, return_interpolated=False, verbose=verbose
                )
            else:
                clean, filt_params = self.noiseDetect_raw(
                    signal, elec_edf, channels_tsv, return_interpolated=False, verbose=verbose
                )

            # Write edfs if requested
            if write_edf_clean:
                # Clean signal
                # convert to list
                clean_list = [clean[i, :] for i in range(clean.shape[0])]
                if os.path.exists(out_edf_path_clean):
                    logging.warning(
                        f"EDF file {out_edf_path_clean} will be overwritten."
                    )
                create_EDF(
                    self.edf_path,
                    out_edf_path_clean,
                    self.processes,
                    chn_labels=chn_labels,
                    signal=clean_list,
                )
                del clean_list
                self.clean_edf = out_edf_path_clean

            if write_edf_int and self.noiseDetect:
                # Interpolated signal
                logging.critical(
                    "Currently, writing the EDF file for the interpolated signal is not supported."
                )
                # convert to list
                # int_list = [interpolated_sig[i,:] for i in range(interpolated_sig.shape[0])]
                # if os.path.exists(out_edf_path_int):
                #     logging.warning(f"EDF file {out_edf_path_int} will be overwritten.")
                # create_EDF(self.edf_path, out_edf_path_int, self.processes, chn_labels = chn_labels, signal = int_list, n_removed = n_removed)
                # del int_list
                # self.inter_edf = out_edf_path_int

            # Write tsv with noise data
            if write_tsv and self.noiseDetect:
                if os.path.exists(out_tsv_path):
                    logging.warning(f" tsv file {out_tsv_path} will be overwritten.")
                df_epochs.to_csv(out_tsv_path, index=False, sep="\t")
            
            # Report
            # First get the table with means of signals
            df_report = pd.DataFrame({
                'Channel': chn_labels,
                'Original mean': np.mean(signal, axis=1),
                'New mean': np.mean(clean, axis=1)
            })
            # Now the report with method for detrending 
            report = {
                'Method used': self.RmTrendMethod,
                'Discarded channels': discarded_labels
            }
            if self.highpass is not None:
                report['Transition band (Hz)'] = self.highpass
            # Finally plot the filter if exists
            if len(filt_params)>0:
                plot_filter(filt_params[0], filt_params[1], os.path.dirname(out_edf_path_clean), self.highpass)
            # Return
            if self.noiseDetect and return_interpolated:
                return clean, interpolated, df_epochs, df_report, report
            elif self.noiseDetect:
                return clean, df_epochs, df_report, report
            # Else, just return the filtered signal
            return clean, df_report, report
                
        except:
            edf_in.close()
            print(traceback.format_exc())
            raise Exception

    def noiseDetect_raw(
        self, raw, electrodes_edf, channels_tsv, return_interpolated=False, verbose=False
    ):
        # electrodes_edf: list with channels present in edf file
        import numpy as np
        import pandas as pd
        import traceback
        import nibabel as nb
        import os

        print(raw.shape)
        if self.noiseDetect and (self.subject == None or self.subjects_dir == None):
            raise Exception(
                "Please indicate a valid subject and directory for the freesurfer outputs."
            )
        # Remove drifts (highpass data) if required
        if self.highpass != None:
            print("Removing trend")
            raw, filt_params = remove_trend(raw, self.RmTrendMethod, self.srate, self.highpass)
            # raw = clean_drifts(raw,self.srate,Transition=self.highpass)
        print(raw.shape)
        # Remove flat-line channels
        # if self.maxFlatlineDuration != None:
        #     print('Removing flatlines')
        #     raw = clean_flatlines(raw,self.srate,maxFlatlineDuration=self.maxFlatlineDuration)
        # print(raw.shape)
        if self.noiseDetect:
            # ---------------Run automatic detection of noise using autoreject-------------------
            print("Running autoreject")
            chn_pos = get_chn_positions(channels_tsv, electrodes_edf, self.tfm)
            # Channels to extract
            keys = list(chn_pos.keys())
            # Number of samples
            N = raw.shape[-1]
            # Create time vector using srate
            t = np.arange(0, N) / self.srate

            # Create sEEG montage
            montage = get_montage(chn_pos, self.subject, self.subjects_dir)
            # Initiate csv epoch file
            cols = ["Epoch #", "Start ID", "End ID"] + keys
            # Create MNE epochs
            mne_epochs, epochs_ids, n_missed = create_mne_epochs(
                raw, keys, self.srate, montage, self.epoch_autoreject
            )
            # Update IDs
            start_IDs = epochs_ids["Start ID"]
            end_IDs = epochs_ids["End ID"]

            # Run autoreject
            epochs_ar, noise_labels = run_autoreject(mne_epochs, verbose=verbose)
            # Create noise df
            # Start-end IDs for each epoch
            IDs_array = np.array([start_IDs, end_IDs]).T
            # Epoch numbering
            epoch_num = np.arange(1, 1 + len(start_IDs))
            noise_array = np.c_[epoch_num, IDs_array, noise_labels]
            df_epochs = pd.DataFrame(data=noise_array, columns=cols)

            # Return interpolated signal only if requested!
            if return_interpolated:
                # Reshape to n_chn x n_time
                clean_sig = epochs_ar.get_data()
                print(clean_sig.shape)
                clean_sig = clean_sig.swapaxes(0, 1).reshape(len(keys), -1)
                # Attach the non-clean part of the signal
                if n_missed != 0:
                    # print(n_missed)
                    # print(signal[:,-n_missed:])
                    sig_missed = raw[:, -n_missed]
                    if sig_missed.ndim == 1:
                        sig_missed = sig_missed.reshape(-1, 1)
                clean_sig = np.hstack([clean_sig, sig_missed])
                return raw, clean_sig, df_epochs, filt_params
            else:
                return raw, df_epochs, filt_params
        else:
            return raw, filt_params
