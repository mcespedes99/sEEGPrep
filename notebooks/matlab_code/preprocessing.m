clear; clc;

filename = '/home/mcesped/scratch/Datasets/bids_ieeg/Subj_space/sub-EPL31LHS0026/ses-V02SE06/ieeg/sub-EPL31LHS0026_ses-V02SE06_task-full_run-01_ieeg.edf';

%% Get the epoch
% define trials
cfg            = [];
cfg.dataset    = filename;
cfg.trialfun   = 'trialfun_annotation';
cfg = ft_definetrial(cfg);  
cfg            = ft_definetrial(cfg);
% read the data
cfg.continuous = 'yes';
% cfg.channel    = {'LPIn2', 'LPIn3'}; %idx=73,74
data           = ft_preprocessing(cfg);
idx = 73;

% Plot 
selected_data = data.trial(1);
selected_data = selected_data{1};
figure;
plot(data.time{1}(1:60*2048), selected_data(idx,1:60*2048))
title('LPIn2 channel plot for the epoch of interest')
xlabel('Time (s)')
ylabel('Amplitude (uV)')

% Save data
signals = data.trial(1);
signals = signals{1};
save('epoching.mat', "signals");

labels = data.label;
save('unipolar_labels_matlab.mat', "labels");

%% Downsample
% the default functionality in ft_resampledata applies a firls
% anti-aliasing filter that has its cutoff at the new Nyquist frequency
cfg = [];
cfg.resamplefs = 200;
data = ft_resampledata(cfg, data);

% Plot 
selected_data = data.trial(1);
selected_data = selected_data{1};
figure;
plot(data.time{1}(1:60*200), selected_data(idx,1:60*200))
title('LPIn2 channel plot - Downsampled')
xlabel('Time (s)')
ylabel('Amplitude (uV)')

% Save data
signals = data.trial(1);
signals = signals{1};
save('Downsampling.mat', "signals");

%% Drift correction using high-pass

cfg = [];
cfg.demean = 'yes';
cfg.baselinewindow = 'all';
cfg.padding = 4*60;
cfg.padtype = 'data';
cfg.hpfilter = 'yes';
cfg.hpfreq = 0.25;
cfg.hpfiltord = 3;
data = ft_preprocessing(cfg, data);

% Plot 
selected_data = data.trial(1);
selected_data = selected_data{1};
figure;
plot(data.time{1}(1:60*200), selected_data(idx,1:60*200))
title('LPIn2 channel plot - Drfit correction')
xlabel('Time (s)')
ylabel('Amplitude (uV)')

% Save data
signals = data.trial(1);
signals = signals{1};
save('Highpass.mat', "signals");


%% Re-reference

depths = {'LOFr*', 'LPCg*', 'LAm*', 'LAHc*', 'LPHc*', 'LTePo*', 'LAIn*', 'LPIn*',...
    'ROFr*', 'RPCg*', 'RAm*', 'RAHc*', 'RPHc*', 'RAIn*', 'RPIn*'};
for d = 1:numel(depths)
cfg = [];
cfg.channel = ft_channelselection(depths{d}, data.label);
cfg.reref = 'yes';
cfg.refchannel = 'all';
cfg.refmethod = 'bipolar';
reref_depths{d} = ft_preprocessing(cfg, data);
end

% Append together
cfg            = [];
data = ft_appenddata(cfg, reref_depths{:});

% Plot 
idx = 56;
selected_data = data.trial(1);
selected_data = selected_data{1};
figure;
plot(data.time{1}(1:60*200), selected_data(idx,1:60*200))
title('LPIn2-3 channel plot - Re-ref')
xlabel('Time (s)')
ylabel('Amplitude (uV)')

% Save data
signals = data.trial(1);
signals = signals{1};
save('Reref.mat', "signals");

labels = data.label;
save('bipolar_labels_matlab.mat', "labels");

%% Run cleanline
data = signals;
signals = [];
signals.data = data;
signals.srate = 200;
lineNoiseIn = struct('lineNoiseMethod', 'clean', ...
        'lineNoiseChannels', 1:size(data,1),...
        'Fs', 200, ...
        'lineFrequencies', [60],...
        'p', 0.01, ...
        'fScanBandWidth', 8, ...
        'taperBandWidth', 2, ...
        'taperWindowSize', 4, ...
        'taperWindowStep', 1, ...
        'tau', 100, ...
        'pad', 2, ...
        'fPassBand', [0 100], ...
        'maximumIterations', 10);
[signal, lineNoiseOut] = cleanLineNoise(signals, lineNoiseIn);
signals = signal.data;
save('cleanline.mat', "signals");