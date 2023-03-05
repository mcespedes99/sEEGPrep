from .utils import get_montage, get_chn_labels, get_orig_data, get_chn_positions, downsampling
from .clean_autoreject import create_mne_epochs, run_autoreject
from .clean_drifts import clean_drifts
from .clean_flatlines import clean_flatlines
from .clean_PLI import removePLI_chns, zapline, cleanline, notch_filt
from .rereference import create_EDF, create_bipolars, extract_location, apply_bipolar_criteria, get_chn_info, create_EDF_from_signal
import logging

class cleanSEEG:
    
    def __init__(self, 
                 edf_path, 
                 chn_csv_path, 
                 subject, 
                 subjects_dir,
                 RmTrendMethod = 'HighPass',
                 cleanPLI = True, 
                 methodPLI = 'Zapline', 
                 lineFreq = 60,
                 bandwidth = 4,
                 n_harmonics = 3,
                 noiseDetect = True,
                 highpass = [0.5, 1], #I set it to [0.5, 1.5] to improve comp cost
                 maxFlatlineDuration = 5, 
                 trsfPath=None, 
                 epoch_length=5, 
                 processes = None):
        import pyedflib
        self.edf_path = edf_path
        self.chn_csv_path = chn_csv_path
        self.subject = subject
        self.subjects_dir = subjects_dir
        self.RmTrendMethod = RmTrendMethod
        self.cleanPLI = cleanPLI
        self.methodPLI = methodPLI
        self.lineFreq = lineFreq
        self.bandwidth = bandwidth
        self.n_harmonics = n_harmonics
        self.noiseDetect = noiseDetect
        self.highpass = highpass # Set to None to shut down
        self.maxFlatlineDuration = maxFlatlineDuration
        self.trsfPath = trsfPath
        self.epoch_length = epoch_length
        self.processes = processes
        # Extra params 
        self.reref_chn_list = []
        self.reref_chn_df = []
        self.inter_edf = None
        self.clean_edf = None
        
        # Find sample rate
        edf_in = pyedflib.EdfReader(edf_path)
        self.srate = edf_in.getSampleFrequencies()[0]/edf_in.datarecord_duration
        edf_in.close()
    
    
    
    def rereference(self,
                    out_edf_path,
                    write_tsv = False,
                    out_tsv_path = None,
                    df_cols = None,
                    use_clean = False):
        # TODO: include use_clean feature!
        import os
        import pyedflib
        import numpy as np
        import pandas as pd
        # Manage some errors:
        if write_tsv and out_tsv_path==None:
            raise Exception('To write the tsv file, please indicate an output path.')
        # TODO: separate rerefering from label map creation
        print('Running rereference')
        # Extract info about electrodes positions
        elec_pos = pd.read_csv(self.chn_csv_path, sep='\t')
        # Create bipolar combinations
        bipolar_channels, bipolar_info_df = create_bipolars(elec_pos, self.processes, df_cols = df_cols)
        # Save values in the class
        self.reref_chn_list = bipolar_channels
        self.reref_chn_df = bipolar_info_df
        # Write tsv
        if write_tsv: #and not os.path.exists(out_tsv_name):
            if os.path.exists(out_tsv_path):
                logging.warning(f"tsv file {out_tsv_path} will be overwritten.")
            bipolar_info_df.to_csv(out_tsv_path, index=False, sep = '\t')
        
        # print(bipolar_channels)
        
        # Create new EDF file
        if os.path.exists(out_edf_path):
            logging.warning(f"edf file {out_edf_path} will be overwritten.")
        create_EDF(self.edf_path, out_edf_path, self.processes, chn_labels = bipolar_channels)
        self.rereference_edf = out_edf_path
    
    def identify_regions(self,
                         aparc_aseg_path,
                         # conf = None,
                         use_reref = True,
                         write_tsv = False,
                         out_tsv_path = None,
                         df_cols = None,
                         use_clean = False,
                         discard_wm_un = False,
                         write_edf = False,
                         out_edf_path = None):
        import os
        # discard_wm_un discards white matter and unknown from edf file (not from csv)
        # Manage a few exceptions
        if write_edf and out_edf_path == None:
            raise Exception('Please indicate a value for out_edf_path or set write_edf to False.')
        
        if use_reref:
            if len(self.reref_chn_list)==0 or len(self.reref_chn_df)==0:
                raise Exception('Please run rereference first or set use_reref to False.')
            chn_info_df = self.reref_chn_df
            df_cols_keys = ['type', 'label', 'x', 'y', 'z', 'group']
            df_cols_vals = chn_info_df.columns.values.tolist()
            df_cols = dict(zip(df_cols_keys, df_cols_vals))
            chn_list = self.reref_chn_list
        # Extract electrodes information if using not rereference results
        else:
            # if conf==None:
            #     raise Exception('Please indicate a proper configuration (unipolar/bipolar).')
            logging.info(f'Extracting channel positions from {self.chn_csv_path}')
            chn_info_df, chn_list, df_cols = get_chn_info(self.chn_csv_path, df_cols = df_cols) #, conf=conf
            # print(chn_list)
        # Create tsv file with information about location of the channels
        df_location = extract_location(aparc_aseg_path, self.trsfPath, chn_info_df, df_cols)
        if write_tsv:
            if os.path.exists(out_tsv_path):
                logging.warning(f" tsv file {out_tsv_path} will be overwritten.")
            df_location.to_csv(out_tsv_path, index=False, sep = '\t')
        
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
                create_EDF(self.edf_path, out_edf_path, self.processes, chn_labels = chn_list)
                self.rereference_edf = out_edf_path
        elif (not discard_wm_un) and write_edf:
            logging.warning(f"EDF file {out_edf_path} will not be written as no updates have been made.")
        return df_location
    
    def downsample(self,
                   target_srate,
                   write_edf = False,
                   out_edf_path = None):
        import pyedflib
        from multiprocessing.pool import Pool
        from multiprocessing import get_context
        from functools import partial
        import numpy as np
        import os
        # Open edf file
        edf_in = pyedflib.EdfReader(self.edf_path)
        # Number of channels to downsample
        n_chns = edf_in.signals_in_file
        edf_in.close()
        # Run through different channels
        channels = np.arange(n_chns)
        # create a process context. Refer to:
        # https://github.com/dask/dask/issues/3759
        ctx = get_context('spawn')
        with Pool(processes=self.processes, context=ctx) as pool:
            data_list, newSrate = zip(*pool.map(partial(downsampling, edf_file=self.edf_path, 
                                                   orig_srate=self.srate, target_srate=target_srate), channels))
        n_samples = len(data_list[0])
        data_dnsampled = np.zeros((n_chns, n_samples))
        for ch in np.arange(len(data_list)):
            data_dnsampled[ch,:] = data_list[ch]
        
        # Write edf if requested
        if write_edf:
            if os.path.exists(out_edf_path):
                logging.warning(f"EDF file {out_edf_path} will be overwritten.")
            create_EDF(self.edf_path, out_edf_path, self.processes, signal = data_list, new_srate = newSrate[0])
            del data_list
            self.downsample_edf = out_edf_path
        return data_dnsampled, newSrate[0]
    
    
    
    def clean_PLI(self, signal):
        if self.methodPLI == 'Cleanline':
            signal = cleanline(signal.T, self.srate, processes = self.processes, bandwidth=self.bandwidth)
        elif self.methodPLI == 'Zapline':
            signal = zapline(signal.T, self.lineFreq/self.srate, self.srate)
        elif self.methodPLI == 'NotchFilter': # add bandwidth param
            signal = notch_filt(signal.T, self.lineFreq, self.srate, n_harmonics = self.n_harmonics)
        elif self.methodPLI == 'PLIremoval':
            signal = removePLI_chns(signal.T, self.srate, 3, [100,0.01,4], [0.1,2,5], 2, processes = self.processes, f_ac=self.lineFreq) #Hardcoded for now
        else:
            raise Exception('PLI method not valid.')
        signal = signal.T
        return signal
    
    # Function to remove trend
    def remove_trend(self, raw):
        import scipy.signal
        import numpy as np
        if self.RmTrendMethod == 'HighPass':
            detsignal = clean_drifts(raw,self.srate,Transition=self.highpass)
        elif self.RmTrendMethod == 'LinearDetrend':
            # linear detrending
            detsignal = scipy.signal.detrend(raw, axis=-1)
        elif self.RmTrendMethod == 'Demean':
            raw_mean = np.mean(raw, axis=-1)
            detsignal = np.subtract(raw, raw_mean.reshape((raw_mean.shape[0],-1)))
        return detsignal
    
    
    def clean_epochs(self,
                     return_interpolated=False, 
                     write_edf_clean = False,
                     out_edf_path_clean = None,
                     write_tsv = False,
                     out_tsv_path = None,
                     write_edf_int = False, # Not working for now
                     out_edf_path_int = None
                    ):
        import pyedflib
        import numpy as np
        import pandas as pd
        import traceback
        import logging
        import os
        # Manage a few exceptions:
        if write_edf_clean and out_edf_path_clean==None:
            raise Exception('EDF file with clean signal cannot be written without and appropiate out path')
        if write_edf_int and out_edf_path_int==None:
            raise Exception('EDF file with interpolated signal cannot be written without and appropiate out path')
        if write_tsv and out_tsv_path==None:
            raise Exception('TSV file with noise information cannot be written without and appropiate out path')
            
        # Begin by getting the position of the electrodes in RAS space
        chn_labels = get_chn_labels(self.chn_csv_path)
        # Extract the labels and timestamps required
        labels, timestamps_epochs = get_orig_data(self.edf_path)
        # Defining start and end time of first epoch
        t_init = timestamps_epochs[0,0]
        t_end = timestamps_epochs[0,1]
        # Open edf file
        edf_in = pyedflib.EdfReader(self.edf_path)
        try:
            # Sample rate
            srate = edf_in.getSampleFrequencies()[0]/edf_in.datarecord_duration
            # Number of samples
            N=edf_in.getNSamples()[0]
            # Create time vector using srate
            t = np.arange(0, N)/self.srate
            ### THIS MUST BE CHANGED TO A MORE STANDARD APPROACH
            # Define length of epochs based on the first one
            t_init_id = np.argmin((np.abs(t-t_init)))
            t_end_id = np.argmin((np.abs(t-t_end)))
            samples_epoch = t_end_id-t_init_id
            print(samples_epoch)
            # Initiate clean signal
            clean = np.array([]).reshape(len(chn_labels),0)
            # Initiate interpolated signal if necessary
            if return_interpolated:
                interpolated_sig = np.array([]).reshape(len(chn_labels),0)
            # Initiate csv epoch file
            cols = ['Epoch #', 'Start ID', 'End ID']+chn_labels
            df_epochs = pd.DataFrame(columns=cols)
            # Last epoch number
            last_epoch = 1
            # Run the algorithm per epoch
            n_epochs = timestamps_epochs.shape[0]
            # Count of removed elements per epoch
            n_removed = [0] 
            # Initialize in 0 as per epoch, the start id will be n_orig[id]-n_removed[id] and the end id: n_orig[id]-n_removed[id+1]
            # Example: if 10 elements were removed in the first epoch and 12 in the 2nd, then the start id of the second should be 
            # the original one minus 10 but the end should be the original minus 22
            for epoch_id in np.arange(n_epochs):
                t_init = timestamps_epochs[epoch_id,0]
                # Find idx for t_init
                t_init_id = np.argmin((np.abs(t-t_init)))
                # Find init for next epoch
                if epoch_id < n_epochs-1:
                    t_init_next = timestamps_epochs[epoch_id+1,0]
                    t_init_next_id = np.argmin((np.abs(t-t_init_next)))
                else:
                    t_init_next_id = N
                # Create signal for that epoch
                signal = np.array([], dtype=np.int64).reshape(0,samples_epoch)
                # Part of the signal not cleaned (zero padded)
                n_not_clean = t_init_next_id - (t_init_id+samples_epoch)
                signal_not_clean = np.array([], dtype=np.int64).reshape(0,n_not_clean)
                # Extract signal per channel
                for chan in chn_labels:
                    id_ch = labels.index(chan)
                    n_extract = t_init_next_id - t_init_id
                    chn_sig = edf_in.readSignal(id_ch, start = t_init_id, n = n_extract)
                    chn_sig_epoch = chn_sig[0:samples_epoch]
                    signal = np.vstack([signal, chn_sig_epoch])
                    signal_not_clean = np.vstack([signal_not_clean, chn_sig[samples_epoch:]])
                edf_in.close()

                # Run cleaning
                # First remove line noise
                print('Removing line noise')
                signal = self.clean_PLI(signal)
                print('PLI removal completed.')
                
                # Second, identify noisy segments and highpass filter the data
                clean_sig, interpolated, tmp_df = self.noiseDetect_raw(signal, t_init_id=t_init_id, epoch_id=last_epoch, return_interpolated=return_interpolated)
                # Attach the non-clean part of the signal
                clean_sig = np.hstack([clean_sig, signal_not_clean])
                interpolated = np.hstack([interpolated, signal_not_clean])
                print(clean_sig.shape)
                
                # Update clean signal
                clean = np.hstack([clean, clean_sig])
                # Update interpolated signal
                interpolated_sig = np.hstack([interpolated_sig, interpolated])
                # Update number of removed elements after autoreject
                n_removed.append(clean.shape[-1]-interpolated_sig.shape[-1])
                # Update output df
                df_epochs = pd.concat([df_epochs, tmp_df])
                # Epochs #s
                last_epoch = last_epoch+len(tmp_df.index)

            # Write edfs if requested
            if write_edf_clean:
                # Clean signal
                # convert to list
                clean_list = [clean[i,:] for i in range(clean.shape[0])]
                if os.path.exists(out_edf_path_clean):
                    logging.warning(f"EDF file {out_edf_path_clean} will be overwritten.")
                create_EDF(self.edf_path, out_edf_path_clean, self.processes, chn_labels = chn_labels, signal = clean_list)
                del clean_list
                self.clean_edf = out_edf_path_clean
            
            if write_edf_int:
                # Interpolated signal
                logging.critical('Currently, writing the EDF file for the interpolated signal is not supported.')
                # convert to list
                # int_list = [interpolated_sig[i,:] for i in range(interpolated_sig.shape[0])]
                # if os.path.exists(out_edf_path_int):
                #     logging.warning(f"EDF file {out_edf_path_int} will be overwritten.")
                # create_EDF(self.edf_path, out_edf_path_int, self.processes, chn_labels = chn_labels, signal = int_list, n_removed = n_removed)
                # del int_list
                # self.inter_edf = out_edf_path_int
            
            # Write tsv with noise data
            if write_tsv:
                if os.path.exists(out_tsv_path):
                    logging.warning(f" tsv file {out_tsv_path} will be overwritten.")
                df_epochs.to_csv(out_tsv_path, index=False, sep = '\t')
            
            # Return 
            if self.noiseDetect and return_interpolated:
                return clean, interpolated_sig, df_epochs
            elif self.noiseDetect:
                return clean, df_epochs
            # Else, just return the filtered signal
            return signal
        except:
            edf_in.close()
            print(traceback.format_exc())
            raise Exception
    
    def noiseDetect_raw(self, raw, t_init_id=0, epoch_id=0, return_interpolated=False):
        import pyedflib
        import numpy as np
        import pandas as pd
        import traceback
        print(raw.shape)
        # Remove drifts (highpass data) if required
        if self.highpass != None:
            print('Removing trend')
            raw = self.remove_trend(raw)
            # raw = clean_drifts(raw,self.srate,Transition=self.highpass)
        print(raw.shape)
        # Remove flat-line channels
        # if self.maxFlatlineDuration != None:
        #     print('Removing flatlines')
        #     raw = clean_flatlines(raw,self.srate,maxFlatlineDuration=self.maxFlatlineDuration)
        # print(raw.shape)
        #---------------Run automatic detection of noise using autoreject-------------------
        print('Running autoreject')
        # Begin by getting the position of the electrodes in RAS space
        chn_pos = get_chn_positions(self.chn_csv_path, self.trsfPath)
        # Channels to extract
        keys = list(chn_pos.keys())
        # Number of samples
        N = raw.shape[-1]
        # Create time vector using srate
        t = np.arange(0, N)/self.srate

        # Create sEEG montage
        montage = get_montage(chn_pos, self.subject, self.subjects_dir)
        # Initiate clean signal
        clean = np.array([]).reshape(len(keys),0)
        # Initiate csv epoch file
        cols = ['Epoch #', 'Start ID', 'End ID']+keys
        # Create MNE epochs
        mne_epochs, epochs_ids, n_missed = create_mne_epochs(raw, keys, self.srate, montage, self.epoch_length)
        # Update IDs
        start_IDs = epochs_ids['Start ID']+t_init_id
        end_IDs = epochs_ids['End ID']+t_init_id

        # Run autoreject
        epochs_ar, noise_labels = run_autoreject(mne_epochs)
        # Create noise df
        # Start-end IDs for each epoch
        IDs_array = np.array([start_IDs,end_IDs]).T
        # Epoch numbering
        epoch_num = np.arange(epoch_id, epoch_id+len(start_IDs))
        noise_array = np.c_[epoch_num, IDs_array, noise_labels]
        df_epochs = pd.DataFrame(data = noise_array, columns = cols)

        # Return interpolated signal only if requested!
        if return_interpolated:
            # Reshape to n_chn x n_time
            clean_sig = epochs_ar.get_data()
            print(clean_sig.shape)
            clean_sig = clean_sig.swapaxes(0,1).reshape(len(keys),-1)
            # Attach the non-clean part of the signal
            if n_missed != 0:
                # print(n_missed)
                # print(signal[:,-n_missed:])
                sig_missed = raw[:,-n_missed]
                if sig_missed.ndim == 1:
                    sig_missed = sig_missed.reshape(-1,1)
            clean_sig = np.hstack([clean_sig, sig_missed])
            return raw, clean_sig, df_epochs
        else:
            return raw, df_epochs

