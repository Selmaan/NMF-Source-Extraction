function [fSub, f_] = removeSourceBaseline_lowpassfilter(f, frameRate)
% Alternative nonparametric baselining function that uses a
% lowpass-filtered version of the fluorescence as the baseline.

if isempty(f) || all(isnan(f))
	warning('Trying to calculate baseline on empty or all NaN vector.')
	fSub = 0;
	f_ = 0;
	return
end

% Specify filter parameters for baselining:
tau = 150; % Filter time constant in seconds.
T = 1/frameRate;
a = T/tau;
nPad = round(1/a)*1;

% Get running percentile as an estimate for the baseline.
f = detrend(f);
n = numel(f);
chunkSize = round(nPad/10);
iStart = 1;
iEnd = chunkSize;
f_ = zeros(1, n);
while iStart < n
    f_(iStart:iEnd) = prctile(f(iStart:min(iEnd, end)), 30);
    iStart = iStart + chunkSize;
    iEnd = iEnd + chunkSize;
end

% Lowpass-filter the running percentile to get the baseline. This can deal
% with baseline flucutuations that are not linear or exponential. Such
% fluctuations seem to be more common in the NMF traces than the manual
% selection.
f_ = padarray(f_, [0 nPad], prctile(f_(1:nPad), 50), 'pre');
f_ = padarray(f_, [0 nPad], prctile(f_(end+1-(1:nPad)), 50), 'post');
f_ = filtfilt(a, [1 a-1], f_);
f_ = f_(nPad + (1:numel(f)));

% Use agressive robustfit to fit the baseline shape to the non-transient
% data points:
warnState = warning('off', 'stats:statrobustfit:IterationLimit');
b = robustfit(f_(:), f(:), 'bisquare', 1);
warning(warnState);
f_ = f_*b(2) + b(1);


fSub = f - f_;