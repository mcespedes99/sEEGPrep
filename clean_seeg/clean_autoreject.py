from .utils import segment_signal

def create_mne_epochs(signal, chn_labels, srate, montage, time_epoch=5):
    import numpy as np
    import mne
    # Divide the signal into small epochs
    signal_epoch, epochs_ids, n_missed = segment_signal(signal, srate, time_epoch)
    # Create information for MNE structure
    info = mne.create_info(ch_names=chn_labels,
                        ch_types=['seeg'] * len(chn_labels),
                        sfreq=srate)
    # Create MNE epoch array 
    mne_epochs = mne.EpochsArray(signal_epoch, info)
    # Set montage
    mne_epochs.set_montage(montage)
    return mne_epochs, epochs_ids, n_missed

def run_autoreject(mne_epoch_array, exclude = []):
    from .autoreject.autoreject import AutoReject, compute_thresholds
    import mne
    import numpy as np
    # Create Autoreject instance
    ar = AutoReject(random_state=42, n_jobs=-1, verbose=True)
    # Run autoreject
    epochs_ar, reject_log = ar.fit_transform(mne_epoch_array, return_log=True)
    # Create 'clean'/'noisy' labels (n_epochs x n_channels)
    noise_labels = np.copy(reject_log.labels).astype('object')
    # Manage possible nan values 
    np.nan_to_num(noise_labels, copy=False)
    # Include bad epochs 
    noise_labels[reject_log.bad_epochs] = 1
    # Define dict for possible values
    noise_map = {0: 'C', 1: 'N', 2: 'N'} #(0.6, 0.6, 0.6, 1.0)
    for key in list(noise_map.keys()):
        noise_labels[noise_labels==key] = noise_map[key]
    return epochs_ar, noise_labels