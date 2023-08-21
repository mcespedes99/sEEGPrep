from .utils import (
    get_montage,
    get_chn_labels,
    get_orig_data,
    get_chn_positions,
    downsampling,
)
from .clean_autoreject import create_mne_epochs, run_autoreject
from .clean_drifts import clean_drifts
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


class cleanSEEG:
    def __init__(
        self,
        edf_path,
        chn_csv_path=None,  # Has to be in MRI RAS space! (not freesurfer)
        RmTrendMethod="HighPass",
        methodPLI="Zapline",
        lineFreq=60,
        bandwidth=4,
        n_harmonics=3,
        noiseDetect=True,
        highpass=[0.5, 1],
        maxFlatlineDuration=5,
        tfm=[],
        epoch_autoreject=5,  # Epoch length for autoreject
        processes=None,
    ):
        import pyedflib

        self.edf_path = edf_path
        self.chn_csv_path = chn_csv_path
        self.RmTrendMethod = RmTrendMethod
        self.methodPLI = methodPLI
        self.lineFreq = lineFreq
        self.bandwidth = bandwidth
        self.n_harmonics = n_harmonics
        self.noiseDetect = noiseDetect
        self.highpass = highpass  # Set to None to shut down
        self.maxFlatlineDuration = maxFlatlineDuration
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
        self, event_label, out_root=None, out_files=None, tmpdir=None, snakemake=False
    ):
        import pyedflib
        import shutil
        import os
        import re
        import snakebids
        import bids
        import numpy as np

        # Find indexes from events
        f = pyedflib.EdfReader(self.edf_path)
        id = [
            value[0]
            for value in enumerate(f.readAnnotations()[2])
            if re.match(event_label, value[1], re.IGNORECASE)
        ]
        # Create df with annotations
        onset_list = f.readAnnotations()[0]
        f.close()
        # Find time stamps where the 'awake trigger' event is happening
        time_stamps = onset_list[id]
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
            entities["suffix"] = "ieeg_tmp.edf"
        else:
            entities["suffix"] = "ieeg.edf"
        del entities["extension"]
        # Add 'task'
        entities["rec"] = "clip"
        # Here call function to create new EDF file
        for index, event_timestamp in enumerate(time_stamps):
            if out_files == None:
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
                create_epoch_EDF(
                    self.edf_path, event_timestamp, out_edf_path, self.processes
                )
            else:
                print("aqui")
                create_epoch_EDF(new_edf, event_timestamp, out_edf_path, self.processes)
        if new_edf != None:
            print("delete")
            os.remove(new_edf)

    # Reference function
    def rereference(
        self,
        out_edf_path,
        write_tsv=False,
        out_tsv_path=None,
        df_cols=None,
        use_clean=False,
    ):
        # TODO: include use_clean feature!
        import os
        import pyedflib
        import numpy as np
        import pandas as pd

        # Manage some errors:
        if write_tsv and out_tsv_path == None:
            raise Exception("To write the tsv file, please indicate an output path.")
        if self.chn_csv_path == None:
            raise Exception(
                "Please indicate the path to a tsv file with channels information."
            )
        # TODO: separate rerefering from label map creation
        print("Running rereference")
        # Extract info about electrodes positions
        elec_pos = pd.read_csv(self.chn_csv_path, sep="\t")
        # Extract channels present in edf file (might not be all)
        edf = pyedflib.EdfReader(self.edf_path)
        elec_edf = edf.getSignalLabels()
        edf.close()
        # Create bipolar combinations
        bipolar_channels, bipolar_info_df = create_bipolars(
            elec_pos, elec_edf, self.processes, df_cols=df_cols
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

    def identify_regions(
        self,
        aparc_aseg_path,
        colortable_file,  # colortable has to be a tsv with at least index, name
        use_reref=True,
        write_tsv=False,
        out_tsv_path=None,
        df_cols=None,
        use_clean=False,
        discard_wm_un=False,
        write_edf=False,
        out_edf_path=None,
        vol_version=False,
        json_out = None
    ):
        import os
        import pyedflib
        import json

        # Manage exception
        if self.chn_csv_path == None:
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
            df_cols_keys = [
                "type",
                "group",
                "label",
                "x_init",
                "x_end",
                "y_init",
                "y_end",
                "z_init",
                "z_end",
            ]
            df_cols_vals = chn_info_df.columns.values.tolist()
            df_cols = dict(zip(df_cols_keys, df_cols_vals))
            chn_list = self.reref_chn_list
        # Extract electrodes information if using not rereference results
        else:
            logging.info(f"Extracting channel positions from {self.chn_csv_path}")
            # First extract list of channels present in edf file
            edf = pyedflib.EdfReader(self.edf_path)
            elec_edf = edf.getSignalLabels()
            edf.close()
            chn_info_df, chn_list, df_cols = get_chn_info(
                self.chn_csv_path, elec_edf, df_cols=df_cols, vol_version=vol_version
            )  # , conf=conf
        # TODO: verify format of df_cols based on "vol_version"
        # Create tsv file with information about location of the channels
        df_location, regions_per_chn = extract_location(
            aparc_aseg_path,
            chn_info_df,
            df_cols,
            self.tfm,
            colortable_file,
            vol_version,
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
                f"EDF file {out_edf_path} will not be written as no updates have been made."
            )
        return df_location

    def downsample(
        self, target_srate, write_edf=False, out_edf_path=None, epoch_edf=False
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
        # Number of channels to downsample
        n_chns = edf_in.signals_in_file
        edf_in.close()
        # Run through different channels
        channels = np.arange(n_chns)
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
        data_dnsampled = np.zeros((n_chns, n_samples))
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
            )
            del data_list
            self.downsample_edf = out_edf_path
        return data_dnsampled, newSrate[0]

    # Function to reject PLI
    def reject_PLI(self, write_edf=False, out_edf_path=None):
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
        # Begin by getting the position of the electrodes in RAS space
        elec_edf = edf_in.getSignalLabels()
        chn_labels = get_chn_labels(self.chn_csv_path, elec_edf)
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

    # Function to remove trend
    def remove_trend(self, raw):
        import scipy.signal
        import numpy as np

        if self.RmTrendMethod == "HighPass":
            detsignal = clean_drifts(raw, self.srate, Transition=self.highpass)
        elif self.RmTrendMethod == "LinearDetrend":
            # linear detrending
            detsignal = scipy.signal.detrend(raw, axis=-1)
        elif self.RmTrendMethod == "Demean":
            raw_mean = np.mean(raw, axis=-1)
            detsignal = np.subtract(raw, raw_mean.reshape((raw_mean.shape[0], -1)))
        return detsignal

    def clean_epochs(
        self,
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
        if self.chn_csv_path == None:
            raise Exception(
                "Please indicate the path to a tsv file with channels information."
            )  # Open edf file
        edf_in = pyedflib.EdfReader(self.edf_path)
        # Extract labels
        elec_edf = edf_in.getSignalLabels()
        # Begin by getting the position of the electrodes in RAS space
        chn_labels = get_chn_labels(self.chn_csv_path, elec_edf)
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
                clean, interpolated, df_epochs = self.noiseDetect_raw(
                    signal, elec_edf, return_interpolated=True, verbose=verbose
                )
                # Update number of removed elements after autoreject
                n_removed.append(clean.shape[-1] - interpolated.shape[-1])
            elif self.noiseDetect:
                clean, df_epochs = self.noiseDetect_raw(
                    signal, elec_edf, return_interpolated=False, verbose=verbose
                )
            else:
                clean = self.noiseDetect_raw(
                    signal, elec_edf, return_interpolated=False, verbose=verbose
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

            # Return
            if self.noiseDetect and return_interpolated:
                return clean, interpolated, df_epochs
            elif self.noiseDetect:
                return clean, df_epochs
            # Else, just return the filtered signal
            return signal
        except:
            edf_in.close()
            print(traceback.format_exc())
            raise Exception

    def noiseDetect_raw(
        self, raw, electrodes_edf, return_interpolated=False, verbose=False
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
            raw = self.remove_trend(raw)
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
            chn_pos = get_chn_positions(self.chn_csv_path, electrodes_edf, self.tfm)
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
                return raw, clean_sig, df_epochs
            else:
                return raw, df_epochs
        else:
            return raw
