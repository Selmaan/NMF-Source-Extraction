function [C,f,A] = pinv_temporal_components(acqObj,nSlice,A,b)

%% Eliminate empty ROIs
emptyROIs = find(sum(A)==0);
if ~isempty(emptyROIs)
    A(:,emptyROIs) = [];
end

%% get activity inference matrix
pA = pinv(full([A,b]));
% pA = pinv_sparse([A,b]);
%% 
imSize = acqObj.correctedMovies.slice(nSlice).channel.size(1,1:2);
nFrames = sum(acqObj.correctedMovies.slice(nSlice).channel.size(:,3));
step = 1000;
frameBatches = 1:step:nFrames;
frameBatches(2,:) = min(frameBatches(1,:) + step - 1, nFrames);
traces = nan(size(A,2)+size(b,2), nFrames);
% traces = cell(1,size(frameBatches,2));
parfor_progress(size(frameBatches,2));
% parfor frameBatch = 1:size(frameBatches,2)
for frameBatch = 1:size(frameBatches,2)
    fInd = frameBatches(1,frameBatch):frameBatches(2,frameBatch);
    dMap = memmapfile(acqObj.indexedMovie.slice(nSlice).channel(1).fileName,...
        'Format', {'int16', [nFrames, prod(imSize)], 'mov'});
    mov = dMap.data.mov;
    tempY = double(mov(fInd,:)');
    tempY = reshape(tempY,imSize(2),imSize(1),length(fInd));
    tempY = permute(tempY,[2 1 3]);
    tempY = reshape(tempY,prod(imSize),length(fInd));
    clear mov
    clear dMap
%     traces{frameBatch} = pA * tempY;
    traces(:,fInd) = pA * tempY;
    parfor_progress;
end
parfor_progress(0);
% traces = [traces{:}];

%%
nSources = size(A,2);
C = double(traces(1:nSources,:));
f = double(traces(nSources+1:end,:));
