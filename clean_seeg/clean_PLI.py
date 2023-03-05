"""
Algorithms to perform power line interference (PLI) removal.
Author: Mauricio Cespedes Tenorio (Western University)
"""
import numpy as np
# removePLI
def removePLI_chns(data, fs, M, B, P, W, processes = None, f_ac=None):
    from multiprocessing.pool import Pool
    from multiprocessing import get_context
    from functools import partial
    # data: n_samples x n_chn
    n_chns = data.shape[-1]
    # Run through different channels
    channels = np.arange(n_chns)
    data_clean = np.zeros(data.shape)
    # create a process context. Refer to:
    # https://github.com/dask/dask/issues/3759
    ctx = get_context('spawn')
    with Pool(processes=processes, context=ctx) as pool:
        data_list = pool.map(partial(removePLI, data = data, fs=fs, M=M, B=B, P=P, W=W, f_ac=f_ac), 
                      channels)
    for ch in np.arange(len(data_list)):
        data_clean[:,ch] = data_list[ch]
    return data_clean

def removePLI(chn, data, fs, M, B, P, W, f_ac = None):
    """
    removePLI Power Line Interference Cancellation 
    Author: Mauricio Cespedes Tenorio <mcespedes99@gmail.com>
    This is the Python version of the code from: https://github.com/mrezak/removePLI 
    which is an implementation of the proposed algorithm in,
    M. R. Keshtkaran and Z. Yang, "A fast, robust algorithm for power line 
    interference cancellation in neural recording," J. Neural Eng., vol. 11,
    no. 2, p. 026017, Apr. 2014.
    http://iopscience.iop.org/1741-2552/11/2/026017
    http://arxiv.org/abs/1402.6862
    Author of original MATLAB code: Mohammad Reza Keshtkaran <keshtkaran.github@gmail.com>

    Usage:
    s = removePLI(x, fs, M, B, P, W, f_ac)
   x, input (contaminated) signal
    s, output (clean) signal
   fs, sample rate in Hz
   M, number of harmonics to remove
   B, contains three elements [B0,Binf,Bst]: 
	- B0, Initial notch bandwidth of the frequency estimator
	- Binf, Asymptotic notch bandwidth of the frequency estimator
	- Bst, Rate of convergence to 95# of the asymptotic bandwidth Binf
   P, contains three elements [P0,np.pinf,Pst]: 
	- P0, Initial settling time of the frequency estimator
	- np.pinf, Asymptotic settling time of the frequency estimator
	- Pst, Rate of convergence to 95# of the asymptotic settling time
	W, Settling time of the amplitude and phase estimator
 	f_ac, Optional argument, the nominal AC frequency if known (50 Hz or 60 HZ)

	EXAMPLE:
		fs = 500
		n = 120*fs #2-min sequence	
		t = 2*np.pi*(1:n)/fs
		fline = 60 + randn #ramdom interference frequency
		s = filter(1,[1,-0.99],100*randn(1,n)) #1/f PSD
		p = 80*sin(fline*t+randn) + 50*sin(2*fline*t+randn)...
		  + 20*sin(3*fline*t+randn) # interference	
		x = s + p
 		sbar = removePLI(x, fs, 3, [100,0.01,4], [0.1,2,5], 3)
 		pwelch(s,[],[],[],fs) title('PSD of the original signal')
 		figure pwelch(x(fs:end),[],[],[],fs) 
       title('PSD of the contaminated signal')
 		figure pwelch(sbar(fs:end),[],[],[],fs) 
     title('PSD after interference cancellation')
    """
    from scipy import signal
    # print(chn)
    x = data[:, chn]
    del data
    x_mean = np.mean(x)
    x = np.subtract(x, x_mean) #removing the mean
    N = len(x)
    s = np.zeros(N)

    # 3dB cutoff bandwidth
    alpha_f = (1-np.arctan(np.pi*B[0]/fs))/(1+np.arctan(np.pi*B[0]/fs))	#initial, \alpha_0
    alpha_inf = (1-np.tan(np.pi*B[1]/fs))/(1+np.tan(np.pi*B[1]/fs)) #asymptotic	
    alpha_st = np.exp(np.log(0.05)/(B[2]*fs+1))	#rate of change

    # frequency estimator's forgetting factors
    lambda_f = np.exp(np.log(0.05)/(P[0]*fs+1))	#initial
    lambda_inf = np.exp(np.log(0.05)/(P[1]*fs+1))	#asymptotic
    lambda_st = np.exp(np.log(0.05)/(P[2]*fs+1))	#rate of change
    # Smoothing parameter (cut-off freq set at 90 Hz)
    gmma = (1-np.tan(0.5*np.pi*min(90,fs/2)/fs))/(1+np.tan(0.5*np.pi*min(90,fs/2)/fs))

    # phase/amplitude estimator forgetting factor
    lambda_a = np.exp(np.log(0.05)/(W*fs+1))
    if type(lambda_a)!=list: lambda_a = lambda_a*np.ones(M)

    # initializing variables
    kappa_f = 0
    kappa_k = np.zeros(M+1)
    D=10
    C=5
    f_n1=0
    f_n2=0

    # -- Alternative initialization:
    #    kappa_f = cos(55*2*np.pi/fs)
    # --

    # initializing the first oscillator
    u_kp = 1*np.ones(M) #u_k
    u_k = 1*np.ones(M) #u'_k

    # initializing the RLS parameters
    r1 = 10*np.ones(M)
    r4 = 10*np.ones(M)
    a = np.zeros(M)
    b = np.zeros(M)


    # IIR Bandpass filtering:
    if(f_ac != None):     # if AC frequency is known
        if type(f_ac)==list and len(f_ac)==2:
            Fc1 = f_ac[0]  # First Cutoff Frequency
            Fc2 = f_ac[1]  # Second Cutoff Frequency               
        else:
            #Custom center frequency of pass band
            Fc1 = f_ac-2  # First Cutoff Frequency
            Fc2 = f_ac+2  # Second Cutoff Frequency
    else:  #if AC frequency is not known
        #Default 40--70 Hz pass band
        Fc1 = 40  # First Cutoff Frequency
        Fc2 = 70  # Second Cutoff Frequency
    

    ordr   = 4   # Order
    frange  = [Fc1, Fc2]
    nyquist = fs/2 #
    fkernB,fkernA = signal.butter(4,np.array(frange)/nyquist,btype='bandpass')
    x_f =  signal.filtfilt(fkernB,fkernA,x)
    x_f = np.concatenate(([0], np.diff(x_f)))		#First Difference
    
    # h  = fdesign.bandpass('N,F3dB1,F3dB2', ordr, Fc1, Fc2, fs)
    # Hd = design(h, 'butter')
    # x_f = filter(Hd,x)	#Bandpass Filtering
    # x_f = [0 diff(x_f)]		#First Difference

    #--------- Start of data processing:
    for n in np.arange(1,N):
        # print(kappa_f)
        # Lattice Filter
        f_n = x_f[n] + kappa_f*(1+alpha_f)*f_n1 - alpha_f*f_n2

        # Frequency Estimation
        C = lambda_f*C+(1-lambda_f)*f_n1*(f_n+f_n2)
        D = lambda_f*D+(1-lambda_f)*2*f_n1**2
        kappa_t=C/D
        if kappa_t <-1: kappa_t=-1
        if kappa_t > 1: kappa_t= 1
        kappa_f = gmma*kappa_f + (1-gmma)*kappa_t
        # print(kappa_f)
        # Updating lattice states
        f_n2=f_n1
        f_n1=f_n  

        # Bandwidth and Forgetting Factor Updates
        alpha_f = alpha_st*alpha_f + (1-alpha_st)*alpha_inf
        lambda_f = lambda_st*lambda_f + (1-lambda_st)*lambda_inf

        # Discrete-Time Oscillators
        kappa_k[1] = kappa_f
        kappa_k[0] = 1

        e=x[n]
        for k in np.arange(0,M): #for each harmonic do:
            # calculating Cos(kw) for k=1,2...
            if(k+2<=M):
                kappa_k[k+2] = 2*kappa_f*kappa_k[k+1] - kappa_k[k];
            kappa = kappa_k[k+1]
            # print(kappa)
            # Oscillator
            tmp = kappa*(u_kp[k]+u_k[k])
            tmp2 =u_kp[k]
            u_kp[k] = tmp - u_k[k]
            u_k[k] = tmp + tmp2

            # Gain Control
            G = 1.5 - (u_kp[k]**2 - (kappa-1)/(kappa+1)*u_k[k]**2)
            if G<=0: G=1
            u_kp[k] = G * u_kp[k]
            u_k[k] = G * u_k[k]

            # Phase/Amplitude Adaptation
            sincmp = a[k]*u_k[k] + b[k]*u_kp[k]
            e = e - sincmp
            # Simplified RLS
            r1[k] = lambda_a[k]*r1[k] + u_k[k]**2
            r4[k] = lambda_a[k]*r4[k] + u_kp[k]**2
            a[k] = a[k] + u_k[k]*e/r1[k]
            b[k] = b[k]  + u_kp[k]*e/r4[k] 
        s[n]=e
        # return None
    # print(kappa_k)
    return s

# Zapline
def zapline(x, fline, srate, nremove=1, p={}, filt=1):
    #[y,yy]=nt_zapline(x,fline,nremove,p,plotflag) - remove power line artifact
    #
    #  y: denoised data
    #  yy: artifact
    #
    #  x: data
    #  fline: line frequency (normalized to sr)
    #  nremove: number of components to remove [default: 1]
    #  p: additional parameters:
    #    p.nfft: size of FFT [default:1024]
    #    p.nkeep: number of components to keep in DSS [default: all]
    #    p.niterations: number of iterations for smoothing filter
    #    p.fig1: figure to use for DSS score [default: 100]
    #    p.fig2: figure to use for results [default: 101]
    #  plotflag: plot
    #
    #Examples:
    #  nt_zapline(x,60/1000) 
    #    apply to x, assuming line frequency=60Hz and sampling rate=1000Hz, plot results
    #  nt_zapline(x,60/1000,4)
    #    same, removing 4 line-dominated components 
    #  p=[];p.nkeep=30; nt_zapline(x,60/1000,4,p);
    #    same, truncating PCs beyond the 30th to avoid overfitting
    #  [y,yy]=nt_zapline(x,60/1000)
    #    return cleaned data in y, noise in yy, don't plot
    #
    import warnings
    import numpy as np
    from . import utils
    if p=={}:
        p['nfft'] = 1024
        p['nkeep'] = []
        p['niterations'] = 1

    # Handling arguments
    if x.size == 0:
        raise Exception("x data cannot be an empty array")
    # Assuming a shape nxm for x:
    if nremove>=x.shape[0]:
        raise Exception("Number of components cannot be larger than lenght of each signal")
    if fline>1/2:
        raise Exception('fline should be less than Nyquist')
    if x.shape[0]<p['nfft']:
        warnings.warn(f'reducing nfft to {str(x.shape[0])}')
        p['nfft']=2*np.floor(x.shape[0]/2)

    if filt==1:
        xx = utils.square_filt(x,1/fline,p['niterations']) # cancels line_frequency and harmonics, light lowpass
    elif filt==2:
        xx = utils.square_notch_filt(x,fline, srate, 3)
    if p['nkeep']==[]: 
        try:
            p['nkeep']=x.shape[1]
        except:
            p['nkeep'] = 1
    # reduce dimensionality to avoid overfitting
    x_rem = x-xx
    c = np.cov(x_rem.T)
    [eigenvalues, eigvecs]= utils.eigen(c)
    # In python, each column of eigvecs represent an eigenvector. So you have to 
    # multiple x * eigvecs. Easy way to say it, a chunk version of eigvecs could
    # be nxk so the only way to multiply it is x*eigvecs
    # This just rotates the data according to the principal components
    xxxx = (x_rem) @ (eigvecs) 

    print('caca')
    # DSS to isolate line components from residual:
    nHarmonics = np.floor((1/2)/fline);
    [c0,c1] = utils.bias_fft(xxxx, fline*np.arange(1,nHarmonics+1), p['nfft']);
    # print('c0')
    # print(c0)
    # print('c1')
    # print(c1)
    print('2')
    todss = utils.dss(c0,c1);
    print('3')
    # This would be the projection of the noise to the main component of the biased
    # noise, which should represent the line noise.
    xxxx= xxxx @ todss[:,0:nremove] # line-dominated components. 
    # return xxxx
    # Denoise
    xxx = utils.denoise_PCA(x-xx,xxxx); # project them out
    del xxxx

    # reconstruct clean signal
    y=xx+xxx
    del xx
    # del xxx
    yy=x-y

    return y

# Cleanline
def cleanline(EEG, srate, processes=None, channels = [], linefreq = [60], p=0.01, bandwidth = 2, taperHalfBW = 1, taperWinSize=4, taperWinStep = 1, tau=100, padding = 2):
    from multiprocessing.pool import Pool
    from multiprocessing import get_context
    from functools import partial
    from .utils import checkTapers, removeLinesMovingWindow
    fPassBand = [0, srate/2]
    maximumIterations = 10
    # Add harmonics:
    if len(linefreq) == 1:
        lf = linefreq[0]
        for i in range(2,4):
            if i*lf < srate/2:
                linefreq.append(i*lf)
        del lf
    print(linefreq)
    # Set channels if not set
    if len(channels)==0:
        sz = EEG.shape
        channels = np.arange(sz[1]) # EEG should be samples x number channels

    # Set up multi-taper parameters
    taperTemplate = [taperHalfBW, taperWinSize, 1]
    Nwin = round(srate*taperWinSize) # number of samples in window
    # Window cannot be longer than time points
    if Nwin>EEG.shape[0]:
        Nwin = EEG.shape[0]
    tapers = checkTapers(taperTemplate, Nwin, srate)

    # Perform the calculation for each channel separately
    channels = np.sort(channels)
    data = np.zeros(EEG.shape)
    # create a process context. Refer to:
    # https://github.com/dask/dask/issues/3759
    ctx = get_context('spawn')
    with Pool(processes=processes, context=ctx) as pool:
        data_list = pool.map(partial(removeLinesMovingWindow, data=EEG, fPassBand=fPassBand, srate=srate,
                               bandwidth=bandwidth, linefreq=linefreq, maximumIterations=maximumIterations,
                               p=p, padding=padding, tapers=tapers, taperWinSize=taperWinSize, 
                               taperWinStep=taperWinStep, tau=tau), 
                      channels)
    for ch in np.arange(len(data_list)):
        data[:,ch] = data_list[ch]
    # for ch in channels: # Change to multiprocessing!
    #     print(ch)
    #     data[:, ch] = removeLinesMovingWindow(EEG[:, ch], fPassBand, srate, bandwidth, linefreq, maximumIterations, p, padding, tapers, taperWinSize, taperWinStep, tau)

    return data
  
# Notch filtering the data
def notch_filt(x, fline, srate, bandwidth = 1, n_harmonics=None, save_fig=False):
    import scipy.signal
    import scipy.fftpack
    import matplotlib.pyplot as plt

    #
    #  y: filtered data
    # 
    #  x: data to filter (n_samples x n_chans)
    #  T: samples, size of window (can be fractionary) 
    #  nIterations: number of iterations of smoothing operation (large --> gaussian kernel)
    #

    # if integ>=size(x,1);
    #     x=repmat(mean(x),[size(x,1),1,1,1]);
    #     return;
    # end

    # Possible harmonics
    if n_harmonics == None:
        n_harmonics = np.floor((srate/2)/fline) # how many times fline fits in nyquist
    frex2notch = np.arange(1, n_harmonics+1)*fline

    ## narrowband filter to remove line noise

    # initialize filtered signal
    datafilt = x

    # loop over frequencies
    for idx,fi in enumerate(np.arange(0,len(frex2notch))):
        print(idx)
        # create filter kernel using firwin (fir1 in MATLAB)
        frange = [frex2notch[fi]-bandwidth/2, frex2notch[fi]+bandwidth/2]
        # Order of the filter
        order  = int( 150*(srate/frange[0]) * (idx+1))
        order  = order + ~order%2

        # filter kernel
        filtkern = scipy.signal.firwin( order,frange,pass_zero=True,fs=srate )
        # recursively apply to data    
        datafilt = scipy.signal.filtfilt(filtkern,1,datafilt, axis=0)
#         if save_fig == True: ## REQUIRES UPDATE!!
#             # save figure of the kernel and its spectral response
        plt.subplot(121)
        plt.plot(filtkern)
        plt.title('Time domain')

        plt.subplot(122)
        plt.plot(np.linspace(0,srate,10000),np.abs(scipy.fftpack.fft(filtkern,10000))**2)
        plt.xlim([min(0,frex2notch[fi]-30), min(frex2notch[fi]+30, srate/2)])
        plt.title('Frequency domain')
        plt.show()

    return datafilt