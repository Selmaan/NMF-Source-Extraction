function [C,f,A] = pinv_temporal_components(acqObj_data,nSlice,A,b)

%% Eliminate empty ROIs
emptyROIs = find(sum(A)==0);
if ~isempty(emptyROIs)
    A(:,emptyROIs) = [];
end

%% get activity inference matrix
pA = pinv(full([A,b]));
% pA = pinv_sparse([A,b]);
%% 
if ~isnumeric(acqObj_data)
    memMap = matfile(acqObj_data.indexedMovie.slice(nSlice).channel(1).memMap);
    imSize = acqObj_data.correctedMovies.slice(nSlice).channel.size(1,1:2);
    nFrames = sum(acqObj_data.correctedMovies.slice(nSlice).channel.size(:,3))/memMap.dsRatio;
    useMemMap = true;
else
    imSize = [size(acqObj_data,1),size(acqObj_data,2)];
    nFrames = size(acqObj_data,3);
    useMemMap = false;
end
step = 1000;
frameBatches = 1:step:nFrames;
frameBatches(2,:) = min(frameBatches(1,:) + step - 1, nFrames);
traces = nan(size(A,2)+size(b,2), nFrames);
% parfor_progress(size(frameBatches,2));
for frameBatch = 1:size(frameBatches,2)
    fInd = frameBatches(1,frameBatch):frameBatches(2,frameBatch);
    if useMemMap
        tempY = memMap.Yr(:,fInd);
    else
        tempY = reshape(acqObj_data(:,:,fInd),prod(imSize),length(fInd));
    end
    traces(:,fInd) = pA * tempY;
%     parfor_progress;
end
% parfor_progress(0);
% traces = [traces{:}];

%%
nSources = size(A,2);
C = double(traces(1:nSources,:));
f = double(traces(nSources+1:end,:));
