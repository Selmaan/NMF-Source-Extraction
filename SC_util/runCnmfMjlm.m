function runCnmfMjlm(acqObj)

% First, need to run motion correction with memmap creation.
addpath(genpath('D:\GitHub\other\NMF-Source-Extraction'))
addpath(genpath('D:\GitHub\other\NMF-Source-Extraction\ca_source_extraction-master\deconvolution'))

strct = load('G:\mcTest\MM104_160730\MM104_160730_main.mat');
acqObj = strct.MM104_160730_main;

%% Load data:
if exist(acqObj.indexedMovie.slice(1).channel(1).memMap,'file')
    memMap = matfile(acqObj.indexedMovie.slice(1).channel(1).memMap);
else
    error('Could not find memMap File'),
end

fprintf('Loading down-sampled data into memory...'),
Y = memMap.Yr;
fprintf('Done! \n'),

%% Extract sources:
acqObj.extractSources(1)
% acqObj.extractSources(1, reshape(Y, 512, 512, size(Y,2)))
clear Y
update_temporal_components_fromTiff(acqObj);

%% Get traces:
% A = load(acqObj.roiInfo.slice.NMF.filename,'A');
% A = A.A;

% Load:
% load('G:\mcTest\MM104_160730\Slice01_patchResults_v170311.mat')
% A = bsxfun(@rdivide,A,sqrt(sum(A.^2)));

% l = clusterSourcesWithCurrentNn(A);

[dF,deconv,denoised,Gs,Lams,A,b,f] = extractTraces_NMF(acqObj);
dF = cell2mat(dF);
A = cell2mat(A);
denoised = cell2mat(denoised);
deconv = cell2mat(deconv);

% save(fullfile(acqObj.defaultDir, 'deconvResults_AR2'), ...
%     'dF', 'deconv', 'denoised', 'Gs', 'Lams', '-v7.3');

figure(2)
clf
hold on
n = 10;
i = mod(i, 20) + 1;
c = lines(20);
plot(dF(1:n, :)' + (0:n-1)*200, 'k')
plot(real(denoised(1:n, :)' + (0:n-1)*200), 'color', c(i, :))
plot(real(1*deconv(1:n, :)' + (0:n-1)*200), 'r--')

%% my baselining:

mmBase = nan(size(dF));
for i = 1:size(dF, 1)
    if ~all(isnan(dF(i, :)))
        mmBase(i, :) = getF_(dF(i, :));
    end
end