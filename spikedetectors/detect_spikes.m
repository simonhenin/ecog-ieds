function [win_ied, h, u, h_baseline] = detect_spikes(data, fs, baseline_interval, filtered_std, unfiltered_std)

win_ied = []; h = []; u = [];
if ~exist('baseline_interval','var')|isempty(baseline_interval),
    baseline_interval = 1:length(data);
end
if ~exist('filtered_std','var')|isempty(filtered_std),
    filtered_std = 3;
end
if ~exist('unfiltered_std','var')|isempty(unfiltered_std),
    unfiltered_std = 3;
end

% [b,a] =butter(4, [40 80]./(fs/2)); % for LF_07
[b,a] =butter(3, [25 80]./(fs/2)); 
dfilt = filtfilt(b,a,data);

% envelope of filtered signal
h = hilbert(dfilt);
h = abs(h);
% h = sqrt(real(h).^2+imag(h).^2);

%envelope of unfiltered signal
u = hilbert(data);
u = abs(u);
% u = sqrt(real(u).^2+imag(u).^2);


% h_baseline = mean(h(baseline_interval));
% u_baseline = mean(u(baseline_interval));
h_baseline = median(h(baseline_interval));
u_baseline = median(u(baseline_interval));

% ied_idx = find(h >= h_baseline*filtered_std & u >= u_baseline*unfiltered_std);
% iidx = find( diff(ied_idx) > 1);
% if isempty(ied_idx)|isempty(iidx),
%     return;
% end
% 
% 
% ied_overlap = round( 1*fs);
% win_ied = [ied_idx(1) ied_idx(iidx(1))];
% data_len = length(data);
% for i=1:length(iidx)-1,
%     
%     if ied_idx(iidx(i)+1) - win_ied(end,2) < ied_overlap,
%         continue;
%     end
%     
%     wwin = ied_idx(iidx(i)+1)-ied_overlap:ied_idx(iidx(i)+1)+ied_overlap;
%     
%     if wwin(1)<=0|wwin(end)>data_len,
%         continue;
%     end
%     
%     try,
%         [~,idx_min] = min(data(wwin));
%     catch,
%         keyboard;
%     end
%     
% 
%     win_ied = [win_ied; ...
%         wwin(idx_min) wwin(idx_min)+ied_overlap];
% %     win_ied = [win_ied; ...
% %         wwin(1) wwin(end)];
%     
% end

tt = (h >= h_baseline*filtered_std);
tt(1) = 0; % make sure it starts off as false
tt = diff(tt);

idx1 = find(tt==1); idx2 = find(tt==-1);
if length(idx1) > length(idx2),
    idx1 = idx1(1:length(idx2));
end
idx = [idx1 idx2];

% remove if unfiltered envelope < 3x
pop = [];
for i=1:size(idx,1),
%     uu = u( idx(i,1):idx(i,2) );
%     if sum(  uu < u_baseline*unfiltered_std ) > 0,
    if u( idx(i,1) ) < u_baseline*unfiltered_std,
        pop = [pop i];
    end
end
idx(pop, :) = [];

% remove ied within 1s of previous event
pop = []; ied_overlap = 1*fs;
for i=2:size(idx,1),
     if idx(i,1) - idx(i-1,1) < ied_overlap,
        pop = [pop i];
    end
end
idx(pop, :) = [];

% centers = [];
% for i=1:size(idx,1),
%     hh = h( idx(i, 1):idx(i,2) );
%     [~,c] = max( hh );
%     centers(i) = idx(i,1)+c-1;
% end

win_ied = idx;


