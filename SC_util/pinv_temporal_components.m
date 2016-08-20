function [C,f,A] = pinv_temporal_components(Y,A,b)

%% Eliminate empty ROIs
emptyROIs = find(sum(A)==0);
if ~isempty(emptyROIs)
    A(:,emptyROIs) = [];
end

%% get activity inference matrix
pA = pinv(full([A,b]));

%% 
sizY = Y.sizY;
nFrames = sizY(end);
d = prod(sizY(1:2));
step = 50;
frameBatches = 1:step:nFrames;
frameBatches(2,:) = min(frameBatches(1,:) + step - 1, nFrames);
% traces = nan(size(pA,1),nFrames);
parfor_progress(size(frameBatches,2));
parfor frameBatch = 1:size(frameBatches,2)
%     fInd = frameBatches(frameBatch)+1 : frameBatches(frameBatch+1);
    fInd = frameBatches(1,frameBatch):frameBatches(2,frameBatch);
    tempY = single(Y.Yr(:,fInd));
    traces{frameBatch} = pA * tempY;
    parfor_progress;
end
parfor_progress(0);
traces = [traces{:}];    
%%
nSources = size(A,2);
C = double(traces(1:nSources,:));
f = double(traces(nSources+1:end,:));
