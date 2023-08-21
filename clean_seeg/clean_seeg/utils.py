"""
Utility functions for iEEGPrep.
Author: Mauricio Cespedes Tenorio
"""
import numpy as np
import pandas as pd
import mne


def downsampling(chn, edf_file, orig_srate, target_srate):
    from scipy.interpolate import griddata
    import scipy.signal
    import pyedflib

    # Based on code from Mike X Cohen course (Complete neural signal processing and analysis: Zero to hero)
    # Get the signal
    edf_in = pyedflib.EdfReader(edf_file)
    signal = edf_in.readSignal(chn)
    edf_in.close()
    if orig_srate < target_srate:
        raise Exception(
            "This function only works to downsample and the given original sample rate is bigger than the target one."
        )
    elif orig_srate == target_srate:
        return signal
    # Get info from given signal
    npnts = signal.shape[-1]
    time = np.arange(0, npnts) / orig_srate

    # Upsample the signal close to a multiple of target_srate
    upsampleFactor = int(np.ceil(orig_srate / target_srate))
    newSrate = upsampleFactor * target_srate
    # need to round in case it's not exact
    newNpnts = np.round(npnts * (newSrate / orig_srate))

    # new time vector after upsampling
    newTime = np.arange(0, newNpnts) / newSrate

    ## continue on to interpolation

    # cut out extra time points: we don't want any data above the original last data
    newTime = newTime[newTime <= time[-1]]

    # the new sampling rate actually implemented
    newSrateActual = np.round(1 / np.mean(np.diff(newTime)), 4)
    # print(newSrateActual)

    # interpolate using griddata
    upsampled_signal = griddata(time, signal, newTime, method="cubic")

    ## Downsample to target_srate
    downsampleFactor = int(np.floor(newSrateActual / target_srate))
    downsampledSrate = newSrateActual / downsampleFactor

    # new time vector after upsampling
    newTv = np.arange(0, newTime[-1], 1 / downsampledSrate)
    newPnts = len(newTv)
    # print(newPnts)

    ### low-pass filter at new Nyquist frequency!
    fkern = scipy.signal.firwin(
        int(14 * downsampledSrate / 2),
        downsampledSrate / 2,
        fs=newSrateActual,
        pass_zero=True,
    )
    fsignal = scipy.signal.filtfilt(fkern, 1, upsampled_signal)

    # now downsample
    signal_dsG = fsignal[:-1:downsampleFactor]
    # print(len(signal_dsG))

    return signal_dsG, downsampledSrate


def get_chn_positions(chn_csv_path, electrodes_edf, tfm_list=None):
    """Creates dictionary with the position of each electrode.
    Parameters
    ----------
    ch_csv_path : str
        Path to csv containing electrodes positions.

    Returns : dictionary with position (x,y,z) of each electrode.
    -------

    """
    elec_pos = pd.read_csv(chn_csv_path, sep="\t")
    chn_pos = {}
    for i in np.arange(len(elec_pos)):
        label = elec_pos.loc[[i], ["label"]].values[0][0]
        if label in electrodes_edf:
            pos = elec_pos.loc[[i], ["x", "y", "z"]].values[0]
            for tfm, inv_bool in tfm_list:
                if type(tfm) == str:
                    tfm = readRegMatrix(tfm)
                if inv_bool:
                    tfm = np.linalg.inv(tfm)
                pos = mne.transforms.apply_trans(tfm, pos)
            pos = (pos / 1000).tolist()
            chn_pos[label] = pos
    return chn_pos


def get_chn_labels(chn_csv_path, electrodes_edf):
    """Gets label for each electrode.
    Parameters
    ----------
    ch_csv_path : str
        Path to csv containing electrodes positions.

    Returns : list with labels.
    -------

    """
    elec_info = pd.read_csv(chn_csv_path, sep="\t")
    labels_csv = elec_info["label"].values.tolist()
    # Filter to get only labels that are also in edf file
    labels = [label for label in labels_csv if label in electrodes_edf]
    return labels


def get_montage(ch_pos, subject, subjects_dir):
    """Get montage given Surface RAS (aka mri coordinates in MNE)
    Recovered from:
    https://mne.discourse.group/t/question-re-seeg-source-localization-using-mne/4411/5
    Parameters
    ----------
    ch_pos : dict
        Dictionary of channel positions. Keys are channel names and values
        are 3D coordinates - array of shape (3,) - in native digitizer space
        in m.
    subject ： str
        The name of subject in FreeSurfer
    subjects_dir : str
        The directory of your FreeSurfer subject directory
    contrast : bool
        Boolean indicating whether it is required to transform the coordinates from
        contrast to non-contrast space.

    Returns : head montage
    -------

    """
    subj_trans = mne.coreg.estimate_head_mri_t(subject, subjects_dir)
    mri_to_head_trans = mne.transforms.invert_transform(subj_trans)
    print("Start transforming mri to head")
    print(mri_to_head_trans)

    montage_mri = mne.channels.make_dig_montage(ch_pos, coord_frame="mri")
    montage = montage_mri.copy()
    montage.add_estimated_fiducials(subject, subjects_dir)
    montage.apply_trans(mri_to_head_trans)
    return montage


def readRegMatrix(trsfPath):
    """Reads transformation matrix.
    Parameters
    ----------
    trsfPath : str
            Path to transform from contrast RAS to non-constrast RAS.

    Returns : Transform from contrast RAS to non-contrast RAS space.
    -------

    """

    with open(trsfPath) as (f):
        return np.loadtxt(f.readlines())


def contrast_to_non_contrast(pnt, tfm):
    """Transforms electrode coordinate from contrast space to non-contrast space.
    Parameters
    ----------
    pnt : ndarray (3,)
        Coordinate of the electrode to transform, in mm and RAS space
        from the contrast MRI.
    tfm : ndarray
        Transform from contrast RAS to non-contrast RAS space.

    Returns : Coordinates in RAS for non-constrast MRI space.
    -------

    """
    mri_ras_mm = mne.transforms.apply_trans(tfm, pnt)
    return mri_ras_mm


def get_orig_data(epoch_path):
    """Get labels and start-end timestamps for each epoch."""
    import pyedflib

    # Read data
    edf_in = pyedflib.EdfReader(epoch_path)
    # Read labels
    labels = edf_in.getSignalLabels()
    # Get annotations from edf
    annot = edf_in.readAnnotations()
    annot = {"Onset": annot[0], "Duration": annot[1], "event": annot[2]}
    annot = pd.DataFrame(annot)

    # Get start and end times
    start_times = annot.Onset.to_numpy()[annot["event"].str.match(r"Epoch #\d starts.")]
    end_times = annot.Onset.to_numpy()[annot["event"].str.match(r"Epoch #\d ends.")]
    # Concatenate two epochs
    timestamps_epochs = np.c_[start_times, end_times]
    edf_in.close()

    return labels, timestamps_epochs


def segment_signal(signal, srate, time_epoch=5):
    n_epoch = int(time_epoch * srate)  # 5 seconds by default
    # Initialize segmented signal
    signal_epoch = np.zeros((int(signal.shape[1] / n_epoch), signal.shape[0], n_epoch))
    n_missed = signal.shape[1] - int(signal.shape[1] / n_epoch)
    id = 0
    start_id = []
    end_id = []
    # Segment the signal
    for epoch_id in np.arange(int(signal.shape[1] / n_epoch)):
        tmp = signal[:, id : id + n_epoch]
        start_id.append(id)
        end_id.append(id + n_epoch)
        signal_epoch[epoch_id, :, :] = tmp
        id += n_epoch
    epochs_ids = {"Start ID": start_id, "End ID": end_id}
    return signal_epoch, epochs_ids, n_missed


# Extra function, not used!
def clean_signal(
    edf_path, chn_csv_path, subject, subjects_dir, trsfPath=None, time_epoch=5
):
    import pyedflib
    import traceback

    # Begin by getting the position of the electrodes in RAS space
    chn_pos = get_chn_positions(chn_csv_path, trsfPath)
    # Extract the labels and timestamps required
    labels, timestamps_epochs = get_orig_data(edf_path)
    # Defining start and end time of first epoch
    t_init = timestamps_epochs[0, 0]
    t_end = timestamps_epochs[0, 1]
    # Open edf file
    edf_in = pyedflib.EdfReader(edf_path)
    try:
        # Channels to extract
        keys = list(chn_pos.keys())
        # Sample rate
        srate = edf_in.getSampleFrequencies()[0] / edf_in.datarecord_duration
        # Number of samples
        N = edf_in.getNSamples()[0]
        # Create time vector using srate
        t = np.arange(0, N) / srate
        # Define length of epochs based on the first one
        t_init_id = np.argmin((np.abs(t - t_init)))
        t_end_id = np.argmin((np.abs(t - t_end)))
        length_epoch = t_end_id - t_init_id
        print(length_epoch)
        # Create sEEG montage
        montage = get_montage(chn_pos, subject, subjects_dir)
        # Initiate clean signal
        clean = np.array([]).reshape(len(keys), 0)
        # Initiate csv epoch file
        cols = ["Epoch #", "Start ID", "End ID"] + keys
        df_epochs = pd.DataFrame(columns=cols)
        # Last epoch number
        last_epoch = 0
        # Run the algorithm per epoch
        n_epochs = timestamps_epochs.shape[0]
        for epoch_id in np.arange(n_epochs):
            t_init = timestamps_epochs[epoch_id, 0]
            # Find idx for t_init
            t_init_id = np.argmin((np.abs(t - t_init)))
            # Find init for next epoch
            if epoch_id < n_epochs - 1:
                t_init_next = timestamps_epochs[epoch_id + 1, 0]
                t_init_next_id = np.argmin((np.abs(t - t_init_next)))
            else:
                t_init_next_id = N
            # Create signal for that epoch
            signal = np.array([], dtype=np.int64).reshape(0, length_epoch)
            n_not_clean = t_init_next_id - (t_init_id + length_epoch)
            signal_not_clean = np.array([], dtype=np.int64).reshape(0, n_not_clean)
            # Extract signal per channel
            for chan in keys:
                id_ch = labels.index(chan)
                n_extract = t_init_next_id - t_init_id
                chn_sig = edf_in.readSignal(id_ch, start=t_init_id, n=n_extract)
                chn_sig_epoch = chn_sig[0:length_epoch]
                signal = np.vstack([signal, chn_sig_epoch])
                signal_not_clean = np.vstack([signal_not_clean, chn_sig[length_epoch:]])
            edf_in.close()
            # Create MNE epochs
            mne_epochs, epochs_ids, n_missed = create_mne_epochs(
                signal, keys, srate, montage, time_epoch
            )
            # Update not clean signal if some segments where missed
            if n_missed != 0:
                # print(n_missed)
                # print(signal[:,-n_missed:])
                sig_missed = signal[:, -n_missed]
                if sig_missed.ndim == 1:
                    sig_missed = sig_missed.reshape(-1, 1)
                signal_not_clean = np.hstack([sig_missed, signal_not_clean])
            # Update IDs
            start_IDs = epochs_ids["Start ID"] + t_init_id
            end_IDs = epochs_ids["End ID"] + t_init_id
            # Epochs #s
            epoch_num = np.arange(last_epoch, last_epoch + len(start_IDs))
            last_epoch = last_epoch + len(start_IDs)
            # Run autoreject
            epochs_ar, noise_labels = run_autoreject(mne_epochs)
            # Create noise df
            IDs_array = np.array([start_IDs, end_IDs]).T
            noise_array = np.c_[epoch_num, IDs_array, noise_labels]
            df_tmp = pd.DataFrame(data=noise_array, columns=cols)
            df_epochs = pd.concat([df_epochs, df_tmp])

            # Reshape to n_chn x n_time
            clean_sig = epochs_ar.get_data()
            print(clean_sig.shape)
            clean_sig = clean_sig.swapaxes(0, 1).reshape(len(keys), -1)
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


""" Zapline tools """


def square_filt(x, T, nIterations=1):
    # y=nt_smooth(x,T,nIterations,nodelayflag) - smooth by convolution with square window
    #
    #  y: smoothed data
    #
    #  x: data to smooth
    #  T: samples, size of window (can be fractionary)
    #  nIterations: number of iterations of smoothing operation (large --> gaussian kernel)
    #
    import scipy.signal

    integ = int(np.floor(T))
    frac = T - integ

    # if integ>=size(x,1);
    #     x=repmat(mean(x),[size(x,1),1,1,1]);
    #     return;
    # end

    # remove onset step
    mn = np.mean(x[0 : integ + 1, :], axis=0)
    x = x - mn

    if nIterations == 1 and frac == 0:
        # faster
        x = np.cumsum(x)
        x[T:, :] = x[T:, :] - x[0:-T, :]
        x = x / T
    else:
        # filter kernel
        B = np.concatenate((np.ones(integ), [frac])) / T
        for k in np.arange(1, nIterations):
            B = np.convolve(B, B)
            print("aqui")
        x = scipy.signal.filtfilt(B, 1, x, axis=0)  # lfilter

    # restore DC
    x = x + mn
    return x


# Possible alternative to square filter
def square_notch_filt(x, fline, srate, nHarmonics, plotting=False):
    # y=nt_smooth(x,T,nIterations,nodelayflag) - smooth by convolution with square window
    #
    #  y: smoothed data
    #
    #  x: data to smooth
    #  T: samples, size of window (can be fractionary)
    #  nIterations: number of iterations of smoothing operation (large --> gaussian kernel)
    #

    # if integ>=size(x,1);
    #     x=repmat(mean(x),[size(x,1),1,1,1]);
    #     return;
    # end
    import numpy as np
    import scipy.signal
    import scipy.fftpack

    fline = fline * srate
    # remove onset step
    mn = np.mean(x, axis=0)
    x = x - mn
    # Apply filter
    lower_trans = 0.1
    upper_trans = 0.1
    norder = 24
    filtorder = norder * np.round(srate / 55) + 1
    for n in np.arange(1, nHarmonics + 1):
        if n > 1:
            filtorder = int(filtorder / (n * 0.6))  # Change 0.6 to a parameter
            if filtorder % 2 == 0:
                filtorder += 1
        f_harmonic = n * fline
        lower_bnd = f_harmonic - 5
        upper_bnd = f_harmonic + 5
        # filter kernel
        filter_shape = [1, 1, 0, 0, 1, 1]
        filter_freqs = [
            0,
            lower_bnd * (1 - lower_trans),
            lower_bnd,
            upper_bnd,
            upper_bnd + upper_bnd * upper_trans,
            srate / 2,
        ]
        filter_kern = scipy.signal.firls(
            filtorder, filter_freqs, filter_shape, fs=srate
        )
        # Apply filter:
        x = scipy.signal.filtfilt(filter_kern, 1, x, axis=0)
        if plotting:
            fig, ax2 = plt.subplots(nrows=1, ncols=1, figsize=(15, 6))
            # Power spectrum
            hz = np.linspace(0, srate / 2, int(np.floor(len(filter_kern) / 2) + 1))
            filterpow = np.abs(scipy.fftpack.fft(filter_kern)) ** 2
            ax2.plot(hz, filterpow[: len(hz)], "ks-")
            plt.plot(filter_freqs, filter_shape, "ro-")
            ax2.set_xlim([0, srate / 2])
            ax2.set_xlabel("Frequency (Hz)")
            ax2.set_ylabel("Filter gain")
            ax2.set_title("Frequency response")
            plt.show()
    # restore DC
    x = x + mn
    return x


def bias_fft(x, freq, nfft):
    # [c0,c1]=nt_bias_fft(x,freq,nfft) - covariance with and w/o filter bias
    #
    # x: data
    # freq: row vector of normalized frequencies to keep (wrt sr)
    # nfft: fft size
    #
    # The filter has zeros at all frequencies except those immediately inferior
    # or superior to values in vector freq.
    #
    # If freq has two rows, keep frequencies between corresponding values on
    # first and second row.
    #
    # NoiseTools
    import scipy.fftpack
    import numpy as np

    if max(freq) > 0.5:
        raise Exception("frequencies should be <= 0.5")
    if nfft > x.shape[0]:
        raise Exception("nfft too large")

    # Here the filter is built in the freq domain, which has a reflection
    # Left half of the filter
    filt = np.zeros(int(np.floor(nfft / 2) + 1))
    for k in np.arange(0, len(freq)):
        idx = int(freq[k] * nfft + 0.5)
        filt[idx] = 1

    filt = np.concatenate((filt, np.flipud(filt[0:-1])))

    ## now for convolution
    n = x.shape[0]
    k = len(filt)
    nConv = n + k - 1
    half_kern = int(np.floor(k / 2))
    filt = scipy.fftpack.fft(scipy.fftpack.ifft(filt), nConv)

    # FFTs
    dataX = scipy.fftpack.fft(x, nConv, axis=0)

    # IFFT
    x_filt = np.multiply(dataX, filt.reshape(len(filt), 1))
    x_filt = np.real(scipy.fftpack.ifft(x_filt, axis=0))
    x_filt = x_filt[half_kern:-half_kern, :]

    c0 = np.cov(x.T)
    c1 = np.cov(x_filt.T)

    # return x_filt
    return c0, c1


def eigen(A):
    # Calculates ordered real part of eigenvalues and eigenvectors
    from numpy.linalg import eig
    import numpy as np

    [eigvals, eigvecs] = eig(A)
    eigvecs = np.real(eigvecs)
    eigvals = np.real(eigvals)
    idx = np.flipud(np.argsort(eigvals))
    eigvals = np.flipud(np.sort(eigvals))
    eigvecs = eigvecs[:, idx]
    return eigvals, eigvecs


def dss(c0, c1):
    # [todss,pwr1,pwr2]=nt_dss0(c0,c1,keep1,keep2) - dss from covariance
    #
    # todss: matrix to convert data to normalized DSS components
    # pwr0: power per component (baseline)
    # pwr1: power per component (biased)
    #
    # c0: baseline covariance
    # c1: biased covariance
    # keep1: number of PCs to retain (default: all)
    # keep2: ignore PCs smaller than keep2 (default: 10.^-9)
    #
    import numpy as np

    if c0.shape != c1.shape:
        raise Exception("C0 and C1 should have same size")
    if c0.shape[0] != c0.shape[1]:
        raise Exception("C0 should be square")

    # Eig vals and vecs from the unbiased covariance
    [eigvals0, eigvecs0] = eigen(c0)
    eigvals0 = np.abs(eigvals0)

    # apply PCA and whitening to the biased covariance
    N = np.diag(np.sqrt(1 / (eigvals0)))
    c2 = np.transpose(N) @ np.transpose(eigvecs0) @ c1 @ eigvecs0 @ N

    # matrix to convert PCA-whitened data to DSS
    [eigvals1, eigvecs1] = eigen(c2)

    # DSS matrix (raw data to normalized DSS)
    todss = eigvecs0 * N * eigvecs1
    N2 = np.diag(np.transpose(todss) @ c0 @ todss)
    todss = todss * np.diag(1 / np.sqrt(N2))  # adjust so that components are normalized
    return todss


def crosscov(x, y):
    import numpy as np

    c = np.transpose(x) @ y
    return c


def denoise_PCA(x, ref):
    mnx = np.mean(x, axis=0)
    x = x - mnx

    mnref = np.mean(ref)
    ref = ref - mnref

    # print(len(x))
    cref = ref.T @ ref
    cref = cref / len(x)

    # The crosscov matrix would be just a way of measuring how each channel of x
    # related to the reference signal
    cxref = crosscov(x, ref)
    cxref = cxref / len(x)

    # regression matrix of x on ref
    # PCA of regressor
    if cref.size > 1:
        print("lolo")
        [eigenvalues, eigvecs] = eigen(cref)
        # cross-covariance between data and regressor PCs
        cxref = np.transpose(cxref)
        r = np.transpose(eigvecs) @ cxref
        eigenvalues = np.reshape(eigenvalues, (len(eigenvalues), 1))
        r = np.multiply(r, 1 / eigenvalues)
        r = eigvecs @ r
    else:
        [eigenvalues, eigvecs] = [1, 1]
        # cross-covariance between data and regressor PCs
        cxref = np.transpose(cxref)
        # print(cxref)
        r = cxref

    # TSPCA
    # Then here r is just how each channel from x (which represents the original noise) variates
    # compare to ref (which is just a mx1 signal that represents the best the PL
    # noise accross channels). So it just means the weights of ref to 'reconstruct'
    # each channel from x. The result of ref*r would be constructing a signal per channel
    # that represents the best x based only in ref, which is the PL noise!!
    z = ref @ r
    # z = z/2
    y = x - z
    mny = np.mean(y, axis=0)
    y = y - mny
    return y


"""Cleanline tools"""


def checkTapers(tapers, N, Fs):
    # Helper function to calculate tapers and, if precalculated tapers are supplied,
    # to check that they (the precalculated tapers) the same length in time as
    # the time series being studied. The length of the time series is specified
    # as the second input argument N. Thus if precalculated tapers have
    # dimensions [N1 K], we require that N1=N.
    # Usage: tapers=dpsschk(tapers,N,Fs)
    # Inputs:
    # tapers        (tapers in the form of:
    #                                   (i) precalculated tapers or,
    #                                   (ii) [NW K] - time-bandwidth product, number of tapers)
    #
    # N             (number of samples)
    # Fs            (sampling frequency - this is required for nomalization of
    #                                     tapers: we need tapers to be such
    #                                     that integral of the square of each taper equals 1
    #                                     dpss computes tapers such that the
    #                                     SUM of squares equals 1 - so we need
    #                                     to multiply the dpss computed tapers
    #                                     by sqrt(Fs) to get the right
    #                                     normalization)
    # Outputs:
    # tapers        (calculated or precalculated tapers)
    # eigs          (eigenvalues)
    import numpy as np
    import scipy.signal as ss

    if len(tapers) == 3:  # Fix taper specification
        # Compute timebandwidth product (half BW * N)
        TW = tapers[0] * tapers[1]
        # Compute number of tapers
        K = int(np.floor(2 * TW - tapers[2]))
        tapers = [TW, K]
    sz = len(tapers)
    if sz == 2:
        tmp, eigs = ss.windows.dpss(N, tapers[0], tapers[1], return_ratios=True)
        # Ideally the tapers shouldn't be multiplied by the sqrt(Fs) here as this is only required for
        # the removeLinesMovingWindow function
        tapers = tmp * np.sqrt(Fs)
        # print(tapers.shape)
    return np.array(tapers).T


def removeLinesMovingWindow(
    channel,
    data,
    fPassBand,
    srate,
    bandwidth,
    linefreq,
    maximumIterations,
    p,
    padding,
    tapers,
    taperWinSize,
    taperWinStep,
    tau,
):
    # Removes significant sine waves from (continuous) data using overlapping windows.
    #
    # Usage:
    #    data = removeLinesMovingWindow(data, lineNoise)
    #
    # Parameters
    #       data        data in [N,1] (a single time column vector)
    #       lineNoise   structure with various parameters set
    #       fPassBand       Frequency band used
    #       srate           Sampling frequency
    #       bandwidth  +/- bandwidth centered on each f0 to scan for significant
    #                       lines (TM)
    #       linefreq Line frequencies to be removed
    #       maximumIterations   Maximum times to iterate removal
    #       p               Significance level cutoff
    #       padding             FFT padding factor
    #       tapers          Precomputed tapers from dpss
    #       taperWinSize Taper sliding window lenght (seconds)
    #       taperWinStep Sliding window step size (seconds)
    #       tau             Window overlap smoothing factor
    #
    # Output:
    #       data           Cleaned up data
    #
    # print(channel)
    data = np.copy(data[:, channel])
    # Window,overlap and frequency information
    if data.ndim == 1:
        data = data[:, np.newaxis]
    N = data.shape[0]
    # print(f'N value: {N}')
    Nwin = int(srate * taperWinSize)  # number of samples in window
    # Window cannot be longer than time points
    if Nwin > N:
        Nwin = N
    Nstep = int(taperWinStep * srate)  # number of samples to step through
    Noverlap = Nwin - Nstep  # number of points in overlap
    x = np.arange(0, Noverlap)
    smooth = 1.0 / (
        1 + np.exp(-tau * (x - Noverlap / 2) / Noverlap)
    )  # sigmoidal function
    smooth = smooth[:, np.newaxis]
    winstart = np.arange(0, N - Nwin, Nstep)
    nw = len(winstart)
    # print(f'value nw: {nw}')
    datafit = np.zeros((winstart[-1] + Nwin, 1))

    fidx = np.zeros(len(linefreq))
    f0 = linefreq
    # [initialSpectrum, f] = calculateSegmentSpectrum(data, lineNoise)
    # initialSpectrum = 10*log10(initialSpectrum)
    # for fk = 1:len(linefreq)
    #     [_, fidx(fk)] = min(abs(f - linefreq(fk)))
    # end

    for iteration in np.arange(maximumIterations):
        f0Mask = bool(np.zeros(len(f0)).tolist)
        # print(iteration)
        for n in np.arange(nw):
            indx = np.arange(winstart[n], (winstart[n] + Nwin))
            datawin = data[indx, :]
            datafitwin, f0Sig = fitSignificantFrequencies(
                datawin, f0, fPassBand, bandwidth, srate, p, padding, tapers
            )
            f0Mask = np.logical_or(f0Mask, f0Sig)
            # datafitwin0 = datafitwin incorrectly placed
            if n > 1:
                datafitwin[0:Noverlap, :] = np.multiply(
                    smooth, datafitwin[0:Noverlap, :]
                ) + np.multiply((1 - smooth), datafitwin0[(Nwin - Noverlap) : Nwin, :])
            # print(f'datafit shape {datafit.shape}')
            datafit[indx, :] = datafitwin
            datafitwin0 = datafitwin  # Moved from above the if statement
        data[0 : len(datafit), :] = data[0 : len(datafit), :] - datafit
    return data.T


def fitSignificantFrequencies(
    data, f0, fPassBand, bandwidth, srate, p, padding, tapers
):
    # Fits significant sine waves to specified peaks in continuous data
    #
    # Usage:
    #   [datafit, f0Significant] = fitSignificantFrequencies(data, f0, lineNoise)
    #
    # Parameters:
    #      data        Single channel -- required
    #      f0          Vector with the line frequencies to be removed
    #      lineNoise   Structure with various parameters set
    #       fPassBand       Frequency band used
    #       bandwidth  +/- bandwidth centered on each f0 to scan for significant
    #                       lines (TM)
    #       srate 	            Sampling frequency
    #       p               Significance level cutoff
    #       padding            FFT padding factor
    #       tapers          Precomputed tapers from dpss
    #
    #  Outputs: datafit, f0Significant, FvalSig, aSig, fSig, sig
    #       datafit          Linear superposition of fitted sine waves
    #       f0Significant    f0 values found to be significant
    #
    N = data.shape[0]

    Fval, A, f, sig = testSignificantFrequencies(
        data, fPassBand, srate, p, padding, tapers
    )
    datafit = np.zeros(N)
    frequencyMask = np.zeros(len(f)).tolist()
    f0Significant = np.zeros(len(f0)).tolist()
    if bandwidth != 0:
        # For each line f0(n), scan f0+-BW/2 for largest significant peak of Fval
        for n in np.arange(len(f0)):
            # Extract scan range around f0 ( f0 +- bandwidth/2 )
            ridx_low = np.argmin(np.abs(f - (f0[n] - bandwidth / 2)), axis=0)
            ridx_high = np.argmin(np.abs(f - (f0[n] + bandwidth / 2)), axis=0)
            Fvalscan = Fval[ridx_low:ridx_high, 0]
            Fvalscan[Fvalscan < sig] = 0
            if np.any(Fvalscan):
                # If there's a significant line, pull the max one
                rmaxidx = np.argmax(Fvalscan, axis=0)
                indx = ridx_low + rmaxidx
                frequencyMask[indx] = 1
                f0Significant[n] = 1
    else:
        # Remove exact lines if significant
        for n in np.arange(len(f0)):
            itemp = np.argmin(np.abs(f - f0[n]))
            frequencyMask[itemp] = Fval[itemp] >= sig
            f0Significant[n] = frequencyMask[itemp]
    # Estimate the contribution of any significant f0 lines
    frequencyMask = list(map(bool, frequencyMask))
    f0Significant = list(map(bool, f0Significant))
    fSig = f[frequencyMask]
    fSig = fSig.reshape(len(fSig), 1)
    aSig = A[frequencyMask]
    aSig = aSig.reshape(len(aSig), 1)
    FvalSig = Fval[frequencyMask]
    if len(fSig) > 0:
        x = np.arange(0, N)
        x = x[:, np.newaxis]
        datafit = np.exp(1j * 2 * np.pi * x @ fSig.T / srate) @ aSig + np.exp(
            -1j * 2 * np.pi * x @ fSig.T / srate
        ) @ np.conj(aSig)
    # print(f'datafit size: {datafit.shape}')
    if datafit.ndim == 1:
        return datafit[:, np.newaxis], f0Significant
    return datafit, f0Significant


def testSignificantFrequencies(data, fPassBand, srate, p, padding, tapers):
    # Computes the F-statistic for sine wave in locally-white noise (continuous data).
    #
    # Usage:
    #     [Fval, A, f, sig ,sd] = testSignificantFrequencies(data,lineNoise)
    #
    # Parameters:
    #      data        Single channel -- required
    #      lineNoise   Structure with various parameters set
    #
    # The lineNoise structure has the following fields set:
    #       fPassBand       Frequency band used
    #       srate 	            Sampling frequency
    #       p               Significance level cutoff
    #       padding            FFT padding factor
    #       tapers          Precomputed tapers from dpss
    #
    #  Outputs:
    #       Fval        F-statistic in frequency x 1 form
    #  	    A		    Line amplitude for X in frequency x 1
    # 	    f		    Frequencies of evaluation
    #       sig         F distribution (1-p)# confidence level
    #
    import scipy.stats

    if len(tapers) == 0:
        raise Exception(
            "testSignificantFrequencies:NoTapers",
            "Must provide a tapers field in the lineNoise parameter structure",
        )
    C = data.shape[1]
    # print(tapers)
    N, K = tapers.shape
    nfft = max(
        int(pow(2, np.ceil(np.log(N) / np.log(2)) + padding)), N
    )  # number of points in fft
    f, findx = getfgrid(srate, nfft, fPassBand)  # frequency grid to be returned
    ## Now compute the taper spectrum
    Kodd = np.arange(1, K, 2)
    Keven = np.arange(0, K, 2)
    J = mtfftc(data, tapers, nfft, srate)  # tapered fft of data - f x K x C
    Jp = J[findx, :, :][
        :, Keven, :
    ]  # drop the even ffts and restrict fft to specified frequency grid - f x K x C
    tapers = np.repeat(
        tapers[:, :, np.newaxis], C, axis=2
    )  # add channel indices to the tapers - t x K x C
    H0 = np.squeeze(
        np.sum(tapers[:, Keven, :], axis=0)
    )  # calculate sum of tapers for even prolates - K x C
    H0 = H0[:, np.newaxis]
    Nf = len(findx)  # Number of frequencies
    H0 = np.repeat(
        H0[:, :, np.newaxis], Nf, axis=2
    )  # Add frequency indices to H0 - K x C x f
    H0 = H0.transpose(
        2, 0, 1
    )  # Permute H0 to get dimensions to match those of Jp - f x K x C

    H0sq = np.sum(
        np.multiply(H0, H0), axis=1
    )  # Sum of squares of H0^2 across taper indices - f x C
    JpH0 = np.sum(
        np.multiply(Jp, H0), axis=1
    )  # sum of the product of Jp and H0 across taper indices - f x C
    A = np.divide(JpH0, H0sq)  # amplitudes for all frequencies and channels
    Kp = Jp.shape[1]  # number of even prolates
    Ap = np.repeat(A[:, :, np.newaxis], Kp, axis=2)  # add the taper index to C
    Ap = Ap.transpose(0, 2, 1)  # permute indices to match those of H0
    Jhat = np.multiply(Ap, H0)  # fitted value for the fft

    num = (K - 1) * np.multiply(np.abs(A) ** 2, H0sq)  # numerator for F-statistic
    den = np.sum(np.abs(Jp - Jhat) ** 2, axis=1) + np.sum(
        np.abs(J[findx, :, :][:, Kodd, :]) ** 2, axis=1
    )  # denominator for F-statistic
    Fval = np.divide(num, den)  # F-statisitic
    sig = scipy.stats.f.ppf(1 - p, 2, 2 * K - 2)  # F-distribution based 1-p# point
    A = A * srate
    return Fval, A, f, sig


def mtfftc(data, tapers, nfft, Fs):
    # Multi-taper fourier transform - continuous data
    #
    # Usage:
    # J=mtfftc(data,tapers,nfft,Fs) - all arguments required
    # Input:
    #       data (in form samples x channels/trials or a single vector)
    #       tapers (precalculated tapers from dpss)
    #       nfft (length of padded data)
    #       Fs   (sampling frequency)
    #
    # Output:
    #       J (fft in form frequency index x taper index x channels/trials)
    from scipy.fft import fft

    NC, C = data.shape  # size of data
    # print(data)
    NK, K = tapers.shape  # size of tapers
    if NK != NC:
        raise Exception("length of tapers is incompatible with length of data")
    tapers = np.repeat(
        tapers[:, :, np.newaxis], C, axis=2
    )  # add channel indices to tapers
    data = np.repeat(data[:, :, np.newaxis], K, axis=2)  # add taper indices to data
    data = data.transpose(
        0, 2, 1
    )  # reshape data to get dimensions to match those of tapers
    data_proj = np.multiply(data, tapers)  # product of data with tapers
    J = fft(data_proj, n=nfft, axis=0) / Fs  # fft of projected data
    return J


def getfgrid(Fs, nfft, fpass):
    # Helper function that gets the frequency grid associated with a given fft based computation
    # Called by spectral estimation routines to generate the frequency axes
    # Usage: [f,findx]=getfgrid(Fs,nfft,fpass)
    # Inputs:
    # Fs        (sampling frequency associated with the data)-required
    # nfft      (number of points in fft)-required
    # fpass     (band of frequencies at which the fft is being calculated [fmin fmax] in Hz)-required
    # Outputs:
    # f         (frequencies)
    # findx     (index of the frequencies in the full frequency grid). e.g.: If
    # Fs=1000, and nfft=1048, an fft calculation generates 512 frequencies
    # between 0 and 500 (i.e. Fs/2) Hz. Now if fpass=[0 100], findx will
    # contain the indices in the frequency grid corresponding to frequencies <
    # 100 Hz. In the case fpass=[0 500], findx=[1 512].
    df = Fs / nfft
    f = np.arange(0, Fs, df)  # all possible frequencies
    f = f[0:nfft]
    if len(fpass) != 1:
        findx = np.argwhere(np.all([f >= fpass[0], f <= fpass[-1]], axis=0)).T
    else:
        findx = np.argmin(np.abs(f - fpass))
    f = f[findx][0]
    return f, np.squeeze(findx.T)
