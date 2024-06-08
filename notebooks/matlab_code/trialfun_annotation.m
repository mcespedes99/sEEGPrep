%% Annotation function
function [trl, event] = trialfun_annotation(cfg)

% read the header, this is needed to determine the sampling rate of EEG channels
hdr = ft_read_header(cfg.dataset);

% read the events, don't detect flanks in a trigger channel but read annotations
event = ft_read_event(cfg.dataset, 'detectflank', []);

% make a selection of the Stimulus annotations
sel = vertcat(event.sample) == 37127083;

% determine the sample numbers of events
smp = [event(sel).sample];

begsample = smp;
endsample = smp+round(4*60*hdr.Fs)-1;
offset    = zeros(size(begsample));

trl = [begsample(:) endsample(:) offset(:)];

% remove trials that overlap with the beginning of the file
sel = trl(:,1)>1;
trl = trl(sel,:);

% remove trials that overlap with the end of the file
sel = trl(:,2)<hdr.nSamples;
trl = trl(sel,:);
end