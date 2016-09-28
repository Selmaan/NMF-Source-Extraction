%% Prepare environment:
addpath(genpath('D:\GitHub\other\NMF-Source-Extraction'));

%% Set up:

acqFile = load('D:\Data\scratch\MM102_160725_main\MM102_160725_main_acq.mat');
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

%% Select Slice and run Patches
extractSourcesNMF(acq,1);