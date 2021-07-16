function [spikeTstart, spikeTend] = spikedetect_wave(x, params, Fe, meanStd)
%[spikeEvents] = SPIKEDETECT_WAVE(x, params, Fe)
%   Interictal epileptic spike detector (Wavelet method)
%   Implementation of method described in :
%   [K. P. Indiradevi, E. Elias, P. S. Sathidevi, S. Dinesh Nayak, and K. Radhakrishnan, 
%   “A multi-level wavelet approach for automatic detection of epileptic 
%   spikes in the electroencephalogram,” 
%   Comput. Biol. Med., vol. 38, pp. 805–816, 2008]
%   Implementation adapted from A.J.E. Geerts "Detection of interictal
%   epileptiform discharges in EEG"
%
%   INPUTS:
%       - x             :   data vector [1*nPnts]
%       - params        :   ? 
%       - Fe            :   Sampling frequency
%       - meanStd       :   Standard deviation averaged across channels
%   OUTPUTS:
%       - spikeEvents   :   Events detected

% Parameters
mergeTimeMax    = 0.5;          % Maximal merging time (s)
minSpikeDuration= 0.05;         % Minimum spike duration (s)
smoothingTime   = 0.05;         % Smoothing time (s)
FeOut           = 256;          % (Hz) - can be 2056 or 1024
thresholdMethod = 'adapted';    % can be static, adapted or both

if FeOut == 256
    wavLevel        = 6;
    wavCoeffs       = [4,5];
elseif FeOut == 1024
    wavLevel        = 7;
    wavCoeffs       = [6,7]; 
else
    error('Wrong value for FeOut');
end

%- Downsample signal to 256 Hz
if Fe<FeOut
    msgbox('Sampling frequency should be 256 Hz or higher');
    return;
end
downsampleFactor = floor(Fe/FeOut);
if downsampleFactor>=2
    x   = decimate(x,downsampleFactor);
end
N           = length(x);
   
%- Wavelet decomposition
[C,L]   = wavedec(x,wavLevel,'db4'); % select one channel

%- Extracting detail coefficients at all scales
[D4,D5] = detcoef(C,L,wavCoeffs);
D4i = interp(D4(3:end-4),2^wavCoeffs(1)).^2;
D5i = interp(D5(3:end-4),2^wavCoeffs(2)).^2;

%- Reconstruct detail coefficients
Y4      = upcoef('d',D4,'db4',wavCoeffs(1),N); 
Y5      = upcoef('d',D5,'db4',wavCoeffs(2),N);

%- Compute thresholds
grwav = 2.2580;
gw4 = grwav/4;
gw5 = grwav/(4*sqrt(2));

thresh4 = (meanStd/gw4)*abs(Y4)*2^4;
thresh5 = (meanStd/gw5)*abs(Y5)*2^5;

thresh4_static = (meanStd*std(D4))*2*wavCoeffs(1);
thresh5_static = (meanStd*std(D5))*2*wavCoeffs(2);


%%  Smooth data
% thresh4 = movingaverage1d(thresh4,round(smoothingTime*FeOut));
% thresh5 = movingaverage1d(thresh5,round(smoothingTime*FeOut));
% D4i     = movingaverage1d(D4i,round(smoothingTime*FeOut));
% D5i     = movingaverage1d(D5i,round(smoothingTime*FeOut));


%% Event selection

switch thresholdMethod
    case 'static'
        sel4ind = D4i>thresh4_static;
        sel5ind = D5i>thresh5_static;
    case 'adapted'
        sel4ind = D4i>thresh4;
        sel5ind = D5i>thresh5;
    case 'both'
        thresh4_static = (meanStd*std(D4))*1.5*wavCoeffs(1);
        thresh5_static = (meanStd*std(D5))*1.5*wavCoeffs(2);
        sel4ind = D4i>thresh4 & D4i>thresh4_static;
        sel5ind = D5i>thresh5 & D5i>thresh5_static;
end

selSpikeInd     = sel4ind & sel5ind;
selSpikeInd(end)= 0;
selSpikeInd     = selSpikeInd(:)';
spikeLimits     = (selSpikeInd - [0,selSpikeInd(1:end-1)]);
spikeIndStarts  = nonzeros((spikeLimits==1) .*(1:length(x)));
spikeIndEnds    = nonzeros((spikeLimits==-1).*(1:length(x)));
%- Convert indice to time 
spikeTstart     = (spikeIndStarts-1)./FeOut;
spikeTend       = (spikeIndEnds-1)./FeOut;

%- Merge close events
[spikeTstart, spikeTend] = mergeevents_time(spikeTstart,spikeTend,mergeTimeMax);

%- Reject short events
eventDuration   = spikeTend-spikeTstart;
eventSel        = eventDuration>minSpikeDuration;
spikeTstart     = spikeTstart(eventSel);
spikeTend       = spikeTend(eventSel);

%% Plot signal, squared detail coefficients and thresholds

% t = 0: (1/FeOut) : (N-1)/FeOut;
% figure;
% ax(1) = subplot(311);
% plot(t,x); axis tight;
% 
% ax(2) = subplot(312);
% plot(t,thresh4,'r'); hold on; axis tight;
% plot(xlim,[thresh4_static,thresh4_static],'k');
% plot(t,D4i); 
% 
% title('Scale 4')
% 
% ax(3) = subplot(313);
% plot(t,thresh5,'r'); hold on; axis tight;
% plot(xlim,[thresh5_static,thresh5_static],'k');
% plot(t,D5i); 
% title('Scale 5')
% 
% linkaxes(ax,'x');
% 
% figure;
% t = 0: (1/256) : (N-1)/256;
% ax(1) = subplot(211);  hold on;
% plot(t,x);
% ax(2) = subplot(212); hold on;
% plot(t,thresh4,'r');
% plot(t,thresh5,'g');
% plot(t,D4i);
% plot(t,D5i,'c');
% legend({'thresh4','thresh5','D4i','D5i'});
% linkaxes(ax,'x');



end

