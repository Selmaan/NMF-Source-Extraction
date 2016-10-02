%% Prepare environment:
addpath(genpath('D:\GitHub\other\NMF-Source-Extraction'));

%% Set up:

acqFile = load('D:\Data\scratch\MM102_160725_main\MM102_160725_main.mat');

acq = acqFile.(char(fieldnames(acqFile)));
acq.defaultDir = 'D:\Data\scratch\MM102_160725_main\';

acq.correctedMovies.slice.channel.fileName = changeDir(...
    acq.correctedMovies.slice.channel.fileName, ...
    acq.defaultDir);

acq.indexedMovie.slice.channel.fileName = changeDir(...
    acq.indexedMovie.slice.channel.fileName, ...
    acq.defaultDir);
acq.indexedMovie.slice.channel.fileName = acq.indexedMovie.slice.channel.fileName{:};

fprintf('%d frames from Acq Obj \n',sum(acq.correctedMovies.slice(end).channel.size(:,3))),
acq.save

% Ensure SI metadata is correct (for some reason, fake metadata is present
% in movies corrected on orchestra...follow up):
fPath = acq.Movies{1};
fPath = strrep(fPath, '/n/data2/hms/neurobio/harvey/matthias', '\\research.files.med.harvard.edu\Neurobio\HarveyLab\Matthias\data\');
fPath = strrep(fPath, '/', '\');
[~, acq.metaDataSI] = tiffRead(fPath, [], 1, 0);

% Create dummy syncInfo:
acq.syncInfo.sliceFrames = sum(acq.correctedMovies.slice.channel.size(:, 3));
acq.syncInfo.validFrameCount = acq.syncInfo.sliceFrames;

%% Load processed data:
load('D:\Data\scratch\MM102_160725_main\Slice01_patchResults_v0913.mat')

%% Select Slice and run Patches
extractSourcesNMF(acq,1);


%% Deconvolve and save clustered sources
cd(acq.defaultDir),
[dF,dF_denoised,dF_deconv,...
    traceBs,traceGs,traceSNs,traceSnScales,A,b,f] = deconv_NMF(acq);
% save('cellData_0921','dF','dF_denoised','dF_deconv',...
%     'traceBs','traceGs','traceSNs','traceSnScales','A','b','f'),

%%
% i = 1;
% i = 502;
i=i+1
% i=i-1

% Get frameRate:
if acqObj.metaDataSI.SI.hFastZ.enable
    frameRate = round(acqObj.metaDataSI.SI.hRoiManager.scanFrameRate...
        /acqObj.metaDataSI.SI.hFastZ.numFramesPerVolume);
else
    frameRate = round(acqObj.metaDataSI.SI.hRoiManager.scanFrameRate);
end

% Specify filter parameters for baselining:
tau = 150; % Filter time constant in seconds.
T = 1/frameRate;
a = T/tau;
nPad = round(1/a)*1;

% Get running percentile as an estimate for the baseline. This should be
% improved and made similar to the custom wfun in getF_. In particular, the
% choice of the percentile will bias the baseline. This should be made
% robust (i.e. ignore transients). Perhaps limit calculation to regions
% where slope of smoothed trace is close to zero?
fluor = dF{:}(i, :);
fluor = detrend(fluor);
n = numel(fluor);
chunkSize = round(nPad/10);
iStart = 1;
iEnd = chunkSize;
prct = zeros(1, n);
while iStart < n
    prct(iStart:iEnd) = prctile(fluor(iStart:min(iEnd, end)), 30);
    iStart = iStart + chunkSize;
    iEnd = iEnd + chunkSize;
end

% Lowpass-filter the running percentile to get the baseline. This can deal
% with baseline flucutuations that are not linear or exponential. Such
% fluctuations seem to be more common in the NMF traces than the manual
% selection.
prct = padarray(prct, [0 nPad], prctile(prct(1:nPad), 50), 'pre');
prct = padarray(prct, [0 nPad], prctile(prct(end+1-(1:nPad)), 50), 'post');
prct = filtfilt(a, [1 a-1], prct);
prct = prct(nPad + (1:numel(fluor)));
