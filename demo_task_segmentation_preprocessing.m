clearvars;

edf_file = '/Users/simonhenin/Data/Liulab/Face-Profession Task/Final/Data/ECOG15_NY733/NY733_FaceProfession_512.EDF';
log_file = '/Users/simonhenin/Data/Liulab/Face-Profession Task/Final/Data/ECOG15_NY733/master_log_scored.mat'

[data, hdr ] = mexSLOAD(edf_file);
fs = hdr.SampleRate;
chs = find(~contains(hdr.Label, {'SG', 'DC', 'EKG'}));
labels = hdr.Label(chs);


%% process the triggers and reconcile with the logfile
load(log_file);
localizer = 1; % was localizer run?

% get number of blocks and num exposures per block
idx_dst = find(strcmp(presentation_order(:, 3), 'DST'));
idx_fix = find(strcmp(presentation_order(:, 3), 'FIX'));
num_blocks = length(idx_dst);
num_per_block = idx_dst-idx_fix-1;

trigger_ch = find(strcmp(hdr.Label, 'DC1'));
figure;
plot(data(:, trigger_ch)); hold on;

% find the triggers
idx_trg = find(data(:, trigger_ch)>4e5); %find all indexes where trigger channel is greater than X


% remove very long "triggers" (e.g. Reconfigure segments)
N=1*fs;
x = diff(idx_trg')==1;
f = find([false,x]~=[x,false]);
g = find(f(2:2:end)-f(1:2:end-1)>=N);
long_t = idx_trg(f(2*g-1));

% %making sure triggers are unique
idx_trg( find(diff(idx_trg)<40)+1) = []; 
[~,idx_long]=ismember(long_t, idx_trg);
idx_trg(idx_long) = [];

% remove testing triggers at the beginning of the run
idx_trg(1:14) = [];
stem(idx_trg, ones(size(idx_trg))*3e6); % makes a plot and plot all the triggers


if localizer,
    % remove localizer indices
    idx_localizer = idx_trg(1:64);
    idx_word_localizer = idx_localizer(5:32);
    idx_face_localizer = idx_localizer(34:63);
    idx_fp = idx_trg(65:end);
else
    idx_fp = idx_trg;
end
h1 = stem(idx_word_localizer, ones(size(idx_word_localizer))*7e6); 
h2 = stem(idx_face_localizer, ones(size(idx_face_localizer))*8e6);
stem(idx_fp, ones(size(idx_fp))*5e6);

 
% loop through each block and get the appropriate trigger starting times
idx_exp = []; idx_cr = [];
pop = [];
for i=1:length(num_per_block),
    if i>1,
        idx_start = sum(num_per_block(1:i-1))*2+(i-1)*3+1;
    else
        idx_start = 1;
    end
   idx_exp{i}  = idx_fp(idx_start+1:idx_start+num_per_block(i));
   idx_cr{i}   = idx_fp(idx_start+num_per_block(i)+2:idx_start+num_per_block(i)+1+num_per_block(i));
   
   pop = [pop idx_start+num_per_block(i)+1]; % get rid of final block trigger

   h3 = stem(idx_exp{i}, ones(size(idx_exp{i}))*5.1e6, 'rx'); % encoding
   h4 = stem(idx_cr{i}, ones(size(idx_cr{i}))*5.2e6, 'ro'); % cue-recall
   %stem(idx_rec{i}, ones(size(idx_rec{i}))*400, 'k');
end
legend([h1 h2 h3 h4], {'Word Localizer', 'Face Localizer', 'Exposure', 'Cued Recall'});
% stem(idx_fp(pop), ones(size(pop))*1e4, 'm');
idx_fp(pop) = [];

master_log = presentation_order;
trial_type = presentation_order(:,3);
CR_score = presentation_order(:,4);


%% filtering and HGP extraction
[b1,a1] = butter(4, [59 61]./hdr.SampleRate*2, 'stop');
[b2,a2] = butter(4, [119 121]./hdr.SampleRate*2, 'stop');
[b3,a3] = butter(4, [179 181]./hdr.SampleRate*2, 'stop');
a                       = poly([roots(a1);roots(a2);roots(a3);]);
b                       = conv(b3, conv(b1, b2));
raw_data = filtfilt(b,a, data(:, chs));


%% extract high gamma power
fband = [70:10:170];
for f=2:length(fband),
    [b, a] = butter(4, [fband(f-1) fband(f)]./fs*2);
    tmp = filtfilt( b ,a, raw_data);
    tmp = abs(hilbert(zscore(tmp)));
    HGP(:, :, f-1) = tmp;
end
HGP = mean(HGP, 3);

%% ied detection
[out, discharges] = spike_detector_hilbert_v23(raw_data, fs, '-k1 3.3 -h 60 -b 15');
spks_out = zeros(size(raw_data));
for j=1:size(out.pos),
    ch = out.chan(j);
    position = round(out.pos(j)*fs);
    spks_out( position, ch) = 1;
end

%% segment the data and put into fieldtrip structures       
data_enc = struct('time', [], 'trial', [], 'trialinfo', [],  'fsample', hdr.SampleRate, 'label', {labels'});
hgp_enc = struct('time', [], 'trial', [], 'trialinfo', [],  'fsample', hdr.SampleRate, 'label', {labels'});
spk_enc = struct('time', [], 'trial', [], 'trialinfo', [],  'fsample', hdr.SampleRate, 'label', {labels'});
data_cue = struct('time', [], 'trial', [], 'trialinfo', [], 'fsample', hdr.SampleRate, 'label', {labels'});
hgp_cue = struct('time', [], 'trial', [], 'trialinfo', [], 'fsample', hdr.SampleRate, 'label', {labels'});
spk_cue = struct('time', [], 'trial', [], 'trialinfo', [], 'fsample', hdr.SampleRate, 'label', {labels'});

prestim = 2; poststim = 6; 
for i=1:length(master_log)
    sample_index  = idx_fp(i); %tells you time index for triggers
    
    %switch trial_type{i}
    switch master_log{i,3}
        
        case {'EXP'}
            interval = sample_index-prestim*fs:sample_index+poststim*fs-1; % take data from the sample 3s after Exposure Trial
            tmp = raw_data( interval, :);
            data_enc.trial{end+1} = tmp';
            data_enc.time{end+1} = (0:length(interval)-1)./fs - prestim;
            data_enc.trialinfo(end+1,1) = CR_score{i};
            
            tmp = HGP( interval, :);
            hgp_enc.trial{end+1} = tmp';
            hgp_enc.time{end+1} = (0:length(interval)-1)./fs - prestim;
            hgp_enc.trialinfo(end+1,1) = CR_score{i};
            
            tmp = spks_out( interval, :);
            spk_enc.trial{end+1} = tmp';
            spk_enc.time{end+1} = (0:length(interval)-1)./fs - prestim;
            spk_enc.trialinfo(end+1,1) = CR_score{i};
            
        case {'CUE'}
            interval = sample_index-prestim*fs:sample_index+poststim*fs-1; % take data from the sample 3s after Cued Recall Trial
            tmp = raw_data( interval, :); 
            data_cue.trial{end+1} = tmp';
            data_cue.time{end+1} = (0:length(interval)-1)./fs - prestim;
            data_cue.trialinfo(end+1,1) = CR_score{i};
            
            tmp = HGP( interval, :);
            hgp_cue.trial{end+1} = tmp';
            hgp_cue.time{end+1} = (0:length(interval)-1)./fs - prestim;
            hgp_cue.trialinfo(end+1,1) = CR_score{i};
            
            tmp = spks_out( interval, :);
            spk_cue.trial{end+1} = tmp';
            spk_cue.time{end+1} = (0:length(interval)-1)./fs - prestim;
            spk_cue.trialinfo(end+1,1) = CR_score{i};
    end
end


%% visualize in fieldtrip
cfg = []; 
cfg.viewmode = 'vertical'; 
ft_databrowser(cfg, data_enc);

%% visual artifact rejection
cfg          = [];
cfg.method   = 'trial';
dummy        = ft_rejectvisual(cfg, data_enc);

%% seisupervised rejection
cfg          = [];
cfg.method   = 'summary';
dummy        = ft_rejectvisual(cfg, data_enc);

%% plot some sample data
cfg = []; 
cfg.trials = (data_cue.trialinfo==1); 
hgp_correct = ft_timelockanalysis(cfg, hgp_cue);
cfg.trials = (data_cue.trialinfo==0); 
hgp_incorrect = ft_timelockanalysis(cfg, hgp_cue);

cfg = [];
cfg.channel       = 'DAMF1';
cfg.baseline      = [-0.25 -0.01];
cfg.xlim = [-0.1 6];
cfg.showlegend    = 'yes';
figure;
ft_singleplotER(cfg, hgp_correct, hgp_incorrect);




