function runSNCpatchNMF(data,saveFile,acqBlocks)

%% Set parameters
% sizY = size(data,'Y');                  % size of data matrix
sizY = data.sizY;
patch_size = [54,54];                   % size of each patch along each dimension (optional, default: [32,32])
overlap = [8,8];                        % amount of overlap in each dimension (optional, default: [4,4])

patches = construct_patches(sizY(1:end-1),patch_size,overlap);
p = 1;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.8;                                  % merging threshold

options = CNMFSetParms(...
    'd1',sizY(1),'d2',sizY(2),...
    'search_method','dilate','se',strel('disk',2,0),...
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'temporal_iter',2,...                       % number of block-coordinate descent steps 
    'ssub',1,...
    'tsub',6,...
    'fudge_factor',0.98,...                     % bias correction for AR coefficients
    'merge_thr',merge_thr,...                    % merging threshold
    'gSig',4,...
    'nB',2, ...
    'medw',[1 1]);

% 'search_method','ellipse','dist',3,...      % search locations when updating spatial components

options.temporal_iterFinal = options.temporal_iter;
%% Get Patches Results

parfor_progress(length(patches));
parfor_progress;
RESULTS = patch_CNMF_SNC(data,patches,p,1);
parfor patchNum = 2:length(patches)
    parfor_progress;
    RESULTS(patchNum) = patch_CNMF_SNC(data,patches,p,patchNum);
end
parfor_progress(0);

%% combine results into one structure
fprintf('Combining results from different patches...');
C = cell2mat({RESULTS(:).C}');
d = prod(sizY(1:2));
A = sparse(d,size(C,1));
P.sn = zeros(sizY(1:2));
P.active_pixels = zeros(sizY(1:2));
IND = zeros(sizY(1:2));
P.b = {};
P.c1 = {};
P.gn = {};
P.neuron_sn = {};
P.psdx = zeros(patches{end}(2),patches{end}(4),size(RESULTS(1).P.psdx,2));
P.pFinal = p;

cnt = 0;
B = sparse(prod(sizY(1:end-1)),length(patches));
MASK = zeros(sizY(1:2));
F = zeros(length(patches),sizY(end));
for i = 1:length(patches)
    for k = 1:size(RESULTS(i).A,2)
            cnt = cnt + 1;
            Atemp = zeros(sizY(1:2));
            Atemp(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4)) = ...
                reshape(RESULTS(i).A(:,k),patches{i}(2)-patches{i}(1)+1,patches{i}(4)-patches{i}(3)+1);            
            A(:,cnt) = sparse(Atemp(:));
    end
    
    b_temp = sparse(sizY(1),sizY(2));
    b_temp(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4)) = reshape(RESULTS(i).b,patches{i}(2)-patches{i}(1)+1,patches{i}(4)-patches{i}(3)+1);  
    B(:,i) = b_temp(:);
    F(i,:) = RESULTS(i).f;
    MASK(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4)) = MASK(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4)) + 1;
    P.sn(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4)) = reshape(RESULTS(i).P.sn,patches{i}(2)-patches{i}(1)+1,patches{i}(4)-patches{i}(3)+1);
    P.active_pixels(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4)) = P.active_pixels(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4)) + ...
        reshape(RESULTS(i).P.active_pixels,patches{i}(2)-patches{i}(1)+1,patches{i}(4)-patches{i}(3)+1);
    IND(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4)) = IND(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4)) + 1;
    P.psdx(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),:) = reshape(RESULTS(i).P.psdx,patches{i}(2)-patches{i}(1)+1,patches{i}(4)-patches{i}(3)+1,[]);
end

A = spdiags(1./MASK(:),0,prod(sizY(1:2)),prod(sizY(1:2)))*A;
B = spdiags(1./MASK(:),0,prod(sizY(1:2)),prod(sizY(1:2)))*B;
% A(A<0) = 0;
% B(B<0) = 0;
% C(C<0) = 0;
% F(F<0) = 0;

fprintf(' done. \n');
clear RESULTS

%% Eliminate Empty Sources and Initialize common background
emptyROIs = find(prod(A==0));
if ~isempty(emptyROIs)
        A(:,emptyROIs) = [];
        C(emptyROIs,:) = [];
end

f = nanmedian(F);
fSplit(1,:) = medfilt1(f,1e3,'truncate');
fSplit(2,:) = f-fSplit(1,:);
fBase = prctile(fSplit(2,:),10);
fSplit(1,:) = fSplit(1,:)+fBase;
fSplit(2,:) = fSplit(2,:)-fBase;
f = fSplit;

%% Merge sources and eliminate very weak sources
fracRetain = 0.9;
nSources = size(A,2);
P.b = cell(nSources,1);
P.c1 = cell(nSources,1);
P.gn = cell(nSources,1);
P.neuron_sn = cell(nSources,1);

Km = 0;
Kn = nSources;

while Km < Kn
    Kn = size(A,2);
    [A,C] = merge_components([],A,[],C,f,P,[],options);
    Km = size(A,2),
end

[A,C] = order_ROIs(A,C);

numRetain = ceil(fracRetain*size(A,2));
A = A(:,1:numRetain);
C = C(1:numRetain,:);
P.b = cell(numRetain,1);
P.c1 = cell(numRetain,1);
P.gn = cell(numRetain,1);
P.neuron_sn = cell(numRetain,1);

%% 1-pass on data to get traces w/out dynamics
% P.p = 0;
% options.temporal_iter = 10;
% [C,f,P] = update_temporal_components(data,A,b,C,f,P,options);


% Remove Source Baseline and neuropil (rough estimate w robust fit)
% This is just to help background + neuropil signals be absorbed into the
% background component by taking (most of) it out of the neural sources
warnState = warning('off', 'stats:statrobustfit:IterationLimit');
for nAcq = 1:size(acqBlocks,1)
    acqInd = acqBlocks(nAcq,1):acqBlocks(nAcq,2);
    fSub = removeSourceBaseline(sum(f(:,acqInd),1));
    parfor nSource = 1:size(C,1)
        warning('off', 'stats:statrobustfit:IterationLimit');
        cSub = removeSourceBaseline(C(nSource,acqInd));
        pilFit = robustfit(fSub,cSub,'bisquare',2);
        C(nSource,acqInd) = ...
            removeSourceBaseline(C(nSource,acqInd) - pilFit(1) - pilFit(2)*fSub);
    end
end
warning(warnState);
 
% Set negative values of C to 0 and rescale C + f
% Rescale because the original lars problem for each pixel minimizes the
% total weight over all sources including background! Don't want to
% penalize the background weight dramatically. This 'shrunk' estimate for
% the background is replaced by OLS at end of algorithm
C(C<0) = 0;
cNorm = sqrt(sum(C.^2,2));
C = bsxfun(@rdivide,C,cNorm);
fNorm = sqrt(sum(f.^2,2));
f = bsxfun(@rdivide,f,fNorm);

fprintf('Updating spatial components... (1)');
warning('off','MATLAB:nargchk:deprecated'),
[A,b,C] = update_spatial_components(data,C,f,A,P,options);
fprintf(' done. \n');

% Re-normalize spatial components
A = bsxfun(@rdivide,A,sqrt(sum(A.^2)));
% scale b as well? maybe not, since it is supposed to be much bigger
% compromise might be scale b so its larger by factor of nSources 
b = bsxfun(@rdivide,b,sqrt(sum(b.^2)))*size(A,2);

fprintf('Updating temporal components... (1)')
[C,f,A] = pinv_temporal_components(data,A,b);
fprintf(' done. \n');

fprintf('Merging overlaping components... (1)')
Km = 0;
Kn = size(A,2);
while Km < Kn
    Kn = size(A,2);
    [A,C,~,~,P] = merge_components([],A,b,C,f,P,[],options);
    Km = size(A,2),
end
fprintf(' done. \n');

% if size(YrA,1) > size(C,1) %If we merged components in line above
%     [C,f,P] = update_temporal_components(data,A,b,C,f,P,options);
% end

%% Save Results
if ~isempty(saveFile)
    save(saveFile,'A','b','C','f','P'),
end

%% update spatial and temporal components

% P.p = 0; % Don't model temporal dynamics on first pass
% options.temporal_iter = 1; %Don't use multiple iterations on first pass
% fprintf('Updating spatial components... (1)');
% warning('off','MATLAB:nargchk:deprecated'),
% [A,b,C] = update_spatial_components(data,C,f,A,P,options);
% fprintf(' done. \n');
% fprintf('Updating temporal components... (1)')
% [C,f,P,S,YrA] = update_temporal_components(data,A,b,C,f,P,options);
% fprintf(' done. \n');
% 
% %% Merge Components
% fprintf('Merging overlaping components... (1)')
% Km = 0;
% Kn = size(A,2);
% 
% while Km < Kn
%     Kn = size(A,2);
%     [A,C,~,~,P,S] = merge_components([],A,b,C,f,P,S,options);
%     Km = size(A,2),
% end
% fprintf(' done. \n');
% 
% %% Update components again and merge
% P.p = P.pFinal; % P.p was set to 0 for first spatio-temporal update
% fprintf('Updating Spatial components (2)... ')
% [A,b,C] = update_spatial_components(data,C,f,A,P,options);
% fprintf(' done. \n');
% fprintf('Updating temporal components (2)... ')
% [C,f,P,S,YrA] = update_temporal_components(data,A,b,C,f,P,options);
% fprintf(' done. \n');
% fprintf('Merging overlaping components... (2)')
% Km = 0;
% Kn = size(A,2);
% while Km < Kn
%     Kn = size(A,2);
%     [A,C,~,~,P,S] = merge_components([],A,b,C,f,P,S,options);
%     Km = size(A,2),
% end
% fprintf(' done. \n');
% %% Final spatio-temporal update
% fprintf('Updating spatial components (3)...');
% options.temporal_iter = options.temporal_iterFinal; %Use multiple iterations for final pass
% [A,b,C] = update_spatial_components(data,C,f,A,P,options);
% fprintf(' done. \n');
% fprintf('Updating temporal components (3)... ')
% [C,f,P,S,YrA] = update_temporal_components(data,A,b,C,f,P,options);
% fprintf(' done. \n');

% %% Save Results
% if ~isempty(saveFile)
%     save(saveFile,'A','b','C','f','P','S','YrA'),
% end