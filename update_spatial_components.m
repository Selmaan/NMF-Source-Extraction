function [A,b,C] = update_spatial_components(binFile,C,f,A_,P,options)

% update spatial footprints and background through Basis Pursuit Denoising
% for each pixel i solve the problem 
%   [A(i,:),b(i)] = argmin sum(A(i,:))
%       subject to || Y(i,:) - A(i,:)*C + b(i)*f || <= sn(i)*sqrt(T);
% for each pixel the search is limited to a few spatial components

% INPUTS:
% Y:    raw data
% C:    temporal components
% f:    temporal background
% A_:   current estimate of spatial footprints (used for determining search locations only) 
% P:    dataset parameters (used for noise values and interpolated entries)

% options    parameter struct (for noise values and other parameters)

% OUTPUTS:
% A:    new estimate of spatial footprints
% b:    new estimate of spatial background
% C:    temporal components (updated only when spatial components are completely removed)

% Written by:
% Eftychios A. Pnevmatikakis, Simons Foundation, 2015

warning('off', 'MATLAB:maxNumCompThreads:Deprecated');
imSize = [options.d1,options.d2];
d = prod(imSize);
T = size(C,2);
K = size(C,1);

if nargin < 6 || isempty(options); options = []; end
if ~isfield(options,'use_parallel'); use_parallel = ~isempty(which('parpool')); else use_parallel = options.use_parallel; end % use parallel toolbox if present
if ~isfield(options,'search_method'); method = []; else method = options.search_method; end     % search method for determining footprint of spatial components
options.sn = P.sn;

IND = determine_search_location(A_(:,1:K),method,options);

Cf = [C;f];
if use_parallel         % solve BPDN problem for each pixel
%     Nthr = max(20*maxNumCompThreads,round(d*T/2^24));
%     Nthr = min(Nthr,round(d/1e3));
    Nthr = 510;
    siz_row = [floor(d/Nthr)*ones(Nthr-mod(d,Nthr),1);(floor(d/Nthr)+1)*ones(mod(d,Nthr),1)];
    indeces = [0;cumsum(siz_row)];
    Yf = cell(Nthr,1);
    A = spalloc(d,size(Cf,1),nnz(IND)+size(Cf,1)*d);
    for nthr = 1:Nthr
        dMap = memmapfile(binFile,'Format', {'int16', [T, d], 'mov'});
        mov = dMap.data.mov;
        [matRow, matCol] = ind2sub(imSize, indeces(nthr)+1:indeces(nthr+1));
        binPixInd = sub2ind([imSize(2), imSize(1)], matCol, matRow);
        Ytemp = double(mov(:, binPixInd))';
        clear mov
        clear dMap
        IND_temp = IND(indeces(nthr)+1:indeces(nthr+1),:);
        Atemp = spalloc(siz_row(nthr),size(Cf,1),nnz(IND_temp));
        Yf{nthr} = Ytemp*f'; 
        %Atemp = Acell{nthr};
        sn_temp = options.sn(indeces(nthr)+1:indeces(nthr+1));
        parfor px = 1:siz_row(nthr)
            fn = ~isnan(Ytemp(px,:));       % identify missing data
            ind = find(IND_temp(px,:));
            if ~isempty(ind);
                ind2 = [ind,K+(1:size(f,1))];
                %[~, ~, a, ~] = lars_regression_noise(Ycell{nthr}(px,fn)', Cf(ind2,fn)', 1, Psnc{nthr}(px)^2*T);
                %[~, ~, a, ~] = lars_regression_noise(Ytemp(px,fn)', Cf(ind2,fn)', 1, Psnc{nthr}(px)^2*T);
                warning('off','MATLAB:nargchk:deprecated'),
                [~, ~, a, ~] = lars_regression_noise(Ytemp(px,fn)', Cf(ind2,fn)', 1, sn_temp(px)^2*T);
                a_sparse = sparse(1,ind2,a');
                Atemp(px,:) = a_sparse';
            end
        end
        if mod(nthr,30) == 0
            fprintf('%2.1f%% of pixels completed \n', indeces(nthr+1)*100/d);
        end
        A(indeces(nthr)+1:indeces(nthr+1),:) = Atemp;
    end
    %A = cell2mat(Acell);
    Yf = cell2mat(Yf);
else
    A = [zeros(d,K),zeros(d,size(f,1))];
    sA = zeros(d1,d2);
    Yf = Y*f';
    for px = 1:d   % estimate spatial components
        fn = ~isnan(Y(px,:));       % identify missing data
        ind = find(IND(px,:));
        if ~isempty(ind);
            ind2 = [ind,K+(1:size(f,1))];
            [~, ~, a, ~] = lars_regression_noise(Y(px,fn)', Cf(ind2,fn)', 1, options.sn(px)^2*T);
            A(px,ind2) = a';
            sA(px) = sum(a);
        end
        if show_sum
            if mod(px,d1) == 0;
               figure(20); imagesc(sA); axis square;  
               title(sprintf('Sum of spatial components (%i out of %i columns done)',round(px/d1),d2)); drawnow;
            end
        end
    end
end

A(isnan(A))=0;
A = sparse(A);
A = threshold_components(A,options);  % post-processing of components

fprintf('Updated spatial components \n');

ff = find(sum(A(:,1:K))==0);           % remove empty components
if ~isempty(ff)
    K = K - length(ff);
    A(:,ff) = [];
    C(ff,:) = [];
end

b = double(max((Yf - A(:,1:K)*double(C(1:K,:)*f'))/(f*f'),0));
A = A(:,1:K);