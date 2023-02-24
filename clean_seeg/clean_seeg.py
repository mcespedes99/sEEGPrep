from .utils import get_montage, get_chn_labels, get_orig_data, get_chn_positions
from .clean_autoreject import create_mne_epochs, run_autoreject

class cleanSEEG:
    
    def __init__(self, edf_path, chn_csv_path, subject, subjects_dir, trsfPath=None, epoch_length=5):
        import pyedflib
        self.edf_path = edf_path
        self.chn_csv_path = chn_csv_path
        self.subject = subject
        self.subjects_dir = subjects_dir
        self.trsfPath = trsfPath
        self.epoch_length = epoch_length
        # Find sample rate
        # Open edf file
        edf_in = pyedflib.EdfReader(edf_path)
        self.srate = edf_in.getSampleFrequencies()[0]/edf_in.datarecord_duration
        edf_in.close()
        
    def clean_epochs(self):
        import pyedflib
        import numpy as np
        import pandas as pd
        import traceback
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
            # Define length of epochs based on the first one
            t_init_id = np.argmin((np.abs(t-t_init)))
            t_end_id = np.argmin((np.abs(t-t_end)))
            samples_epoch = t_end_id-t_init_id
            print(samples_epoch)
            # Initiate clean signal
            clean = np.array([]).reshape(len(chn_labels),0)
            # Initiate csv epoch file
            cols = ['Epoch #', 'Start ID', 'End ID']+chn_labels
            df_epochs = pd.DataFrame(columns=cols)
            # Last epoch number
            last_epoch = 1
            # Run the algorithm per epoch
            n_epochs = timestamps_epochs.shape[0]
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
                clean_sig, tmp_df = self.clean_raw(signal, t_init_id=t_init_id, epoch_id=last_epoch)

                # Update output df
                df_epochs = pd.concat([df_epochs, tmp_df])
                # Epochs #s
                last_epoch = last_epoch+len(tmp_df.index)
                # Attach the non-clean part of the signal
                clean_sig = np.hstack([clean_sig, signal_not_clean])
                print(clean_sig.shape)
                # Update clean signal
                clean = np.hstack([clean, clean_sig])
                print(clean.shape)
            return clean, df_epochs
        except:
            edf_in.close()
            print(traceback.format_exc())
            raise Exception
    
    def clean_raw(self, raw, t_init_id=0, epoch_id=0):
        import pyedflib
        import numpy as np
        import pandas as pd
        import traceback
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

        return clean_sig, df_epochs 
        

def clean_raw(raw, srate, chn_csv_path, subject, subjects_dir, t_init_id=0, epoch_id=0, trsfPath=None, time_epoch=5):
    import pyedflib
    import numpy as np
    import pandas as pd
    import traceback
    # Begin by getting the position of the electrodes in RAS space
    chn_pos = get_chn_positions(chn_csv_path, trsfPath)
    # Channels to extract
    keys = list(chn_pos.keys())
    # Number of samples
    N = raw.shape[-1]
    # Create time vector using srate
    t = np.arange(0, N)/srate

    # Create sEEG montage
    montage = get_montage(chn_pos, subject, subjects_dir)
    # Initiate clean signal
    clean = np.array([]).reshape(len(keys),0)
    # Initiate csv epoch file
    cols = ['Epoch #', 'Start ID', 'End ID']+keys
    # Create MNE epochs
    mne_epochs, epochs_ids, n_missed = create_mne_epochs(raw, keys, srate, montage, time_epoch)
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

    return clean_sig, df_epochs

def clean_epochs(edf_path, chn_csv_path, subject, subjects_dir, trsfPath=None, time_epoch=5):
    import pyedflib
    import numpy as np
    import pandas as pd
    import traceback
    # Begin by getting the position of the electrodes in RAS space
    chn_labels = get_chn_labels(chn_csv_path)
    # Extract the labels and timestamps required
    labels, timestamps_epochs = get_orig_data(edf_path)
    # Defining start and end time of first epoch
    t_init = timestamps_epochs[0,0]
    t_end = timestamps_epochs[0,1]
    # Open edf file
    edf_in = pyedflib.EdfReader(edf_path)
    try:
        # Sample rate
        srate = edf_in.getSampleFrequencies()[0]/edf_in.datarecord_duration
        # Number of samples
        N=edf_in.getNSamples()[0]
        # Create time vector using srate
        t = np.arange(0, N)/srate
        # Define length of epochs based on the first one
        t_init_id = np.argmin((np.abs(t-t_init)))
        t_end_id = np.argmin((np.abs(t-t_end)))
        length_epoch = t_end_id-t_init_id
        print(length_epoch)
        # Initiate clean signal
        clean = np.array([]).reshape(len(chn_labels),0)
        # Initiate csv epoch file
        cols = ['Epoch #', 'Start ID', 'End ID']+chn_labels
        df_epochs = pd.DataFrame(columns=cols)
        # Last epoch number
        last_epoch = 1
        # Run the algorithm per epoch
        n_epochs = timestamps_epochs.shape[0]
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
            signal = np.array([], dtype=np.int64).reshape(0,length_epoch)
            # Part of the signal not cleaned (zero padded)
            n_not_clean = t_init_next_id - (t_init_id+length_epoch)
            signal_not_clean = np.array([], dtype=np.int64).reshape(0,n_not_clean)
            # Extract signal per channel
            for chan in chn_labels:
                id_ch = labels.index(chan)
                n_extract = t_init_next_id - t_init_id
                chn_sig = edf_in.readSignal(id_ch, start = t_init_id, n = n_extract)
                chn_sig_epoch = chn_sig[0:length_epoch]
                signal = np.vstack([signal, chn_sig_epoch])
                signal_not_clean = np.vstack([signal_not_clean, chn_sig[length_epoch:]])
            edf_in.close()
            
            # Run cleaning
            clean_sig, tmp_df = clean_raw(signal, srate, chn_csv_path, subject, subjects_dir, t_init_id=t_init_id, epoch_id=last_epoch, trsfPath=None, time_epoch=5)
            
            # Update output df
            df_epochs = pd.concat([df_epochs, tmp_df])
            # Epochs #s
            last_epoch = last_epoch+len(tmp_df.index)
            # Attach the non-clean part of the signal
            clean_sig = np.hstack([clean_sig, signal_not_clean])
            print(clean_sig.shape)
            # Update clean signal
            clean = np.hstack([clean, clean_sig])
            print(clean.shape)
        return clean, df_epochs
    except:
        edf_in.close()
        print(traceback.format_exc())
        return None