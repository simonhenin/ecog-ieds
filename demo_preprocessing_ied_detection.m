
%% add some folders
addpath('biosig4matlab');
addpath('spikedetectors');
addpath('eegplotfunc');

%% load the data
[data, hdr] = mexSLOAD('/Users/simonhenin/OneDrive - NYU Langone Health/X~ X_ea75ab26-67ac-485f-87b5-e73759b5ca06.EDF');
chs = 1:find( strcmp(hdr.Label, 'SG1'))-1;
data = data(:, chs);
labels = hdr.Label(chs);
% add trailing 0
labels = regexprep(labels, '([A-Z]{1,})(\d{1,1})($)', '$10$2');


%% filter the 60-Hz data
[b1,a1] = butter(4, [59 61]./hdr.SampleRate*2, 'stop');
[b2,a2] = butter(4, [119 121]./hdr.SampleRate*2, 'stop');
[b3,a3] = butter(4, [179 181]./hdr.SampleRate*2, 'stop');
a                       = poly([roots(a1);roots(a2);roots(a3);]);
b                       = conv(b3, conv(b1, b2));
data = filtfilt(b,a, data);


%% plot the data
eegplot(data', 'srate', fs, 'spacing', 150, 'labels', labels);

%% bipolar montage?
[montage, data] = getbipolarmontage(labels, [], data);
labels = montage.labelnew;


%%
fs = hdr.SampleRate;
[out, discharges] = spike_detector_hilbert_v23(data, fs, '-k1 3.65 -h 60 -b 15');


%% look at single-channel detections
spks_out = zeros(size(data));
for j=1:size(out.pos),
    %         ch = chs( out.chan(j) );    % relative to analyzed channels
    
    %         % if spike within 50ms of common average spike, then continue
    %         if ~isempty( find( abs(outca.pos- out.pos(j) ) < 0.05) ),
    %             continue;
    %         end
    
    ch = out.chan(j);
    position = round(out.pos(j)*fs);
    spks_out( position, ch) = 1;
end
eegplot(data', 'srate', fs, 'data2', spks_out'.*500, 'spacing', 150, 'labels', labels);

%% multichannel detections
spks_mv = zeros(size(data));
events = [];
for j=1:size(discharges.MP,1),
    idx = find( discharges.MV(j,:) );
    if length(idx)>=2, % only keep multichannel with 2 or more unambigious detections
        ch = idx;
        position = round( discharges.MP(j, idx(1))*fs );
        spks_mv(position, ch) = 1;
    end
end
eegplot(data', 'srate', fs, 'data2', spks_mv'.*500, 'spacing', 150, 'labels', labels);


%% Line-length detector (Kleen Lab)
spks_ll = zeros(size(data));
for ch=1:size(data,2),
    ets=LLspikedetector(data(:, ch),fs);
    for j=1:length(ets),
        % mark the entire spike
        spks_ll( ets(j,1):ets(j, 2), ch ) = 1;
    end
end
eegplot(data', 'srate', fs, 'data2', spks_ll'.*500, 'spacing', 150, 'labels', labels);

%% run it on all channels
[ets,ech]=LLspikedetector(data',fs);
spks_ll = zeros(size(data));
for j=1:length(ets),
    % mark the entire spike
    spks_ll( ets(j,1):ets(j, 2), ech(j, :) ) = 1;
end
eegplot(data', 'srate', fs, 'data2', spks_ll'.*500, 'spacing', 150, 'labels', labels);




