import numpy as np
from sklearn.metrics import mean_squared_error
from scipy.fft import fft, fftfreq
from scipy.signal import welch

def calculate_snr(clean_signal, filtered_signal):
    # Power of the clean signal
    power_clean_signal = np.mean(clean_signal ** 2)
    
    # Power of the noise (difference between the clean and filtered signals)
    noise = clean_signal - filtered_signal
    power_noise = np.mean(noise ** 2)
    
    # Calculate SNR
    snr = 10 * np.log10(power_clean_signal / power_noise)
    return snr

# Get bandpower per freq band
def bandpower(Pxx, f, fmin, fmax):
    ind_min = np.argmax(f > fmin) - 1
    ind_max = np.argmax(f > fmax) - 1
    return np.trapz(Pxx[ind_min: ind_max], f[ind_min: ind_max], axis=-1)

def compute_metrics(clean_signal, filtered_signal, noisy_signal, srate, linefreq=60):
    ## TIME DOMAIN METRICS
    # SNR
    snr = calculate_snr(clean_signal, filtered_signal)
    
    # RMSE
    rmse_time = np.sqrt(mean_squared_error(clean_signal, filtered_signal))
    
    # Correlation
    r = np.corrcoef(clean_signal, filtered_signal)[0,1]

    ## FREQ DOMAIN METRICS
    length_segment = 3.0
    f, clean_psd = welch(clean_signal, fs=srate, nperseg=int(length_segment*srate))
    f, filtered_psd = welch(filtered_signal, fs=srate, nperseg=int(length_segment*srate))
    f, noisy_psd = welch(noisy_signal, fs=srate, nperseg=int(length_segment*srate))
    
    # Spectral error
    rmse_spectral = np.sqrt(mean_squared_error(clean_psd, filtered_psd))

    ## POWER METRICS
    # Attenuation in dB
    # First we need to look for the band of interest
    attenuation = []
    nharmonics = (srate//2)//linefreq
    for n in range(1,nharmonics+1):
        fmin = linefreq*n-5
        fmax = linefreq*n+5
        filtered_power = bandpower(filtered_psd, f, fmin, fmax)
        noisy_power = bandpower(noisy_psd, f, fmin, fmax)
        attenuation.append(20*np.log10((filtered_power/noisy_power)))
    attenuation = np.mean(attenuation)

    return {
        'SNR (dB)': snr,
        'RMSE - Time Domain ($\mu$V)': rmse_time,
        'Correlation coefficient': r,
        'Spectral error ($\mu$V$^2$/Hz)': rmse_spectral,
        'Attenuation (dB)': attenuation
    }