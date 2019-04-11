function [c, s, b, g, lam, active_set] = harvey_constrained_oasisAR1(y, g, sn, optimize_b,...
    optimize_g, decimate, maxIter, tau_range)

% This is a lightly modified version of the constrained-AR1 deconvolution code
% It has been altered to include 'GetSn' and 'oasisAR1' functions, to
% requires an initial decay parameter (rather than use autocorrelation), to
% add zeros to the spike output argument if it is shorter than the denoised
% fluorescence output (happens sometimes if algo thinks no spikes
% happened).

%% Infer the most likely discretized spike train underlying an AR(1) fluorescence trace
% Solves the sparse non-negative deconvolution problem
%  min 1/2|c-y|^2 + lam |s|_1 subject to s_t = c_t-g c_{t-1} >=s_min or =0

%% inputs:
%   y:  T*1 vector, One dimensional array containing the fluorescence intensities
%withone entry per time-bin.
%   g:  scalar, Parameter of the AR(1) process that models the fluorescence ...
%impulse response.
%   sn:  scalar, standard deviation of the noise distribution
%   optimize_b: bool, optimize baseline if True
%   optimize_g: bool, optimize decay if True
%   decimate: int, decimation factor for estimating hyper-parameters faster
%       on decimated data
%   maxIter:  int, maximum number of iterations
%   tau_range: [tau_min, tau_max], the range of tau

%% outputs
%   c: T*1 vector, the inferred denoised fluorescence signal at each time-bin.
%   s: T*1 vector, discetized deconvolved neural activity (spikes)
%   b: scalar, fluorescence baseline
%   g: scalar, parameter of the AR(1) process
%   lam: scalar, sparsity penalty parameter
%   active_set: npool x 4 matrix, active sets

%% Authors: Pengcheng Zhou, Carnegie Mellon University, 2016
% ported from the Python implementation from Johannes Friedrich

%% References
% Friedrich J et.al., NIPS 2016, Fast Active Set Method for Online Spike Inference from Calcium Imaging


%% input arguments
y = reshape(y, [], 1);
T = length(y);

if ~exist('sn', 'var') || isempty(sn)
    sn = GetSn(y, [1/3, 1/2]);
end
if ~exist('lam', 'var') || isempty(lam);   lam = 0; end
if ~exist('optimize_b', 'var') || isempty(optimize_b)
    optimize_b = true;
end
if ~exist('optimize_g', 'var') || isempty(optimize_g)
    optimize_g = true;
end
if ~exist('decimate', 'var') || isempty(decimate)
    decimate = 1;
else
    decimate = max(1, round(decimate));
end
if ~exist('maxIter', 'var') || isempty(maxIter)
    maxIter = 25;
end
if ~exist('tau_range', 'var') || isempty(tau_range)
   g_range = [0, 1];
else
   g_range = exp(-1./tau_range); 
   g = min(max(g, g_range(1)), g_range(2)); 
end
thresh = sn * sn * T;
lam = 0;

% change parameters due to downsampling
if decimate>1
    decimate = 1;  %#ok<NASGU>
    disp('to be done');
    %     fluo = y;
    %     y = resample(y, 1, decimate);
    %     g = g^decimate;
    %     thresh = thresh / decimate / decimate;
    %     T = length(y);
end
g_converged = false;

%% optimize parameters
tol = 1e-4;
% flag_lam = true;
if ~optimize_b   %% don't optimize the baseline b
    %% initialization
    b = 0;
    [solution, spks, active_set] = oasisAR1(y, g, lam);
    
    %% iteratively update parameters lambda & g
    for miter=1:maxIter
        % update g
        if and(optimize_g, ~g_converged)
            g0 = g;
            [solution, active_set, g, spks] = update_g(y, active_set,lam, g_range);
            if abs(g-g0)/g0 < 1e-3  % g is converged
                g_converged = true;
            end
        end
        res = y - solution;
        RSS = res' * res;
        len_active_set = size(active_set, 1);
        if RSS>thresh  % constrained form has been found, stop
            break;
        else
            % update lam
            update_phi;
            lam = lam + dphi;
        end
    end
else
    %% initialization
    b = quantile(y, 0.15);
    [solution, spks, active_set] = oasisAR1(y-b, g, lam);
    update_lam_b;
    
    %% optimize the baseline b and dependends on the optimized g too
    g_converged = false;
    for miter=1:maxIter
        res = y - solution - b;
        RSS = res' * res;
        len_active_set = size(active_set,1);
        
        if or(abs(RSS-thresh) < tol, sum(solution)<1e-9)
            break;
        else
            %% update b & lamba
            update_phi();
            update_lam_b();
            % update b and g
            % update b and g
            if and(optimize_g, ~g_converged)
                g0 = g;
                [solution, active_set, g, spks] = update_g(y-b, active_set,lam, g_range);
                if abs(g-g0)/g0 < 1e-4
                    g_converged = true;
                end
            end
            
        end
    end
    
end
c = solution;
s = spks;
if length(s)<length(c)
    fprintf('Deconv appended %d entries w/ zero value \n', length(c)-length(s)),
    s(end+1:length(c)) = 0;
end

%% nested functions
    function update_phi()  % estimate dphi to match the thresholded RSS
        zeta = zeros(size(solution));
        maxl = max(active_set(:, 4));
        h = g.^(0:maxl);
        for ii=1:len_active_set
            ti = active_set(ii, 3);
            li = active_set(ii, 4);
            idx = 0:(li-1);
            if ii<len_active_set
                zeta(ti+idx) = (1-g^li)/ active_set(ii,2) * h(1:li);
            else
                zeta(ti+idx) = 1/active_set(ii,2) * h(1:li);
            end
        end
        
        if optimize_b
            zeta = zeta - mean(zeta);
            tmp_res = res - mean(res);
            aa = zeta' * zeta;
            bb = tmp_res' * zeta;
            cc = tmp_res'*tmp_res - thresh;
            dphi = (-bb + sqrt(bb^2-aa*cc)) / aa;
        else
            aa = zeta'*zeta;
            bb = res'*zeta;
            cc = RSS-thresh;
            dphi = (-bb + sqrt(bb^2-aa*cc)) / aa;
        end
        if imag(dphi)>1e-9
            flag_phi = false;
            return;
        else
            flag_phi = true;
        end
        active_set(:,1) = active_set(:,1) - dphi*(1-g.^active_set(:,4));
        [solution, spks, active_set] = oasisAR1([], g, lam, [], active_set);
    end

    function update_lam_b() % estimate lambda  & b
        db = mean(y-solution) - b;
        b = b + db;
        dlam = -db/(1-g);
        
        lam = max(0, lam + dlam);
        % correct the last pool
        active_set(end,1) = active_set(end,1) - lam*g^(active_set(end,4));
        ti = active_set(end,3); li = active_set(end,4); idx = 0:(li-1);
        solution(ti+idx) = max(0, active_set(end,1)/active_set(end,2)) * (g.^idx);
    end

end

%update the AR coefficient: g
function [c, active_set, g, s] = update_g(y, active_set, lam, g_range)
%% inputs:
%   y:  T*1 vector, One dimensional array containing the fluorescence intensities
%withone entry per time-bin.
%   active_set: npools*4 matrix, previous active sets.
%   lam:  scalar, curret value of sparsity penalty parameter lambda.
%   g_range: [g_min, g_max] 

%% outputs
%   c: T*1 vector
%   s: T*1 vector, spike train
%   active_set: npool x 4 matrix, active sets
%   g: scalar

%% Authors: Pengcheng Zhou, Carnegie Mellon University, 2016
% ported from the Python implementation from Johannes Friedrich

%% References
% Friedrich J et.al., NIPS 2016, Fast Active Set Method for Online Spike Inference from Calcium Imaging

%% initialization
len_active_set = size(active_set, 1);  %number of active sets
y = reshape(y,[],1);    % fluorescence data
maxl = max(active_set(:, 4));   % maximum ISI
c = zeros(size(y));     % the optimal denoised trace
if ~exist('g_range', 'var') || isempty(g_range)
   g_range = [0, 1];  
end
%% find the optimal g and get the warm started active_set
g = fminbnd(@rss_g, g_range(1), g_range(2));
yp = y - lam*(1-g);
for m=1:len_active_set
    tmp_h = exp(log(g)*(0:maxl)');   % response kernel
    tmp_hh = cumsum(h.*h);        % hh(k) = h(1:k)'*h(1:k)
    li = active_set(m, 4);
    ti = active_set(m, 3);
    idx = ti:(ti+li-1);
    active_set(m,1) = (yp(idx))'*tmp_h(1:li);
    active_set(m,2) = tmp_hh(li);
end
[c,s,active_set] = oasisAR1(y, g, lam, [], active_set);

%% nested functions
    function rss = rss_g(g)
        h = exp(log(g)*(0:maxl)');   % response kernel
        hh = cumsum(h.*h);        % hh(k) = h(1:k)'*h(1:k)
        yp = y - lam*(1-g);     % include the penalty term
        for ii=1:len_active_set
            li = active_set(ii, 4);
            ti = active_set(ii, 3);
            idx = ti:(ti+li-1);
            tmp_v = max(yp(idx)' * h(1:li) / hh(li), 0);
            c(idx) = tmp_v*h(1:li);
        end
        res = y-c;
        rss = res'*res;     % residual sum of squares
    end
end

%% Auxilliary Functions (copied from separate files in original repo)
function sn = GetSn(Y, range_ff, method)
%% Estimate noise standard deviation

%% inputs:
%   Y: N X T matrix, fluorescence trace
%   range_ff : 1 x 2 vector, nonnegative, max value <= 0.5, range of frequency (x Nyquist rate) over which the spectrum is averaged
%   method: string, method of averaging: Mean, median, exponentiated mean of logvalues (default)

%% outputs:
%   sn: scalar, std of the noise

%% Authors: Pengcheng Zhou, Carnegie Mellon University, 2016
% adapted from the MATLAB implemention by Eftychios Pnevmatikakis and the
% Python implementation from Johannes Friedrich

%% References
% Pnevmatikakis E. et.al., Neuron 2016, Simultaneous Denoising, Deconvolution, and Demixing of Calcium Imaging Data

%% input arguments
if ~exist('range_ff', 'var') || isempty(range_ff)
    range_ff = [.25, .5];
end
if ~exist('method', 'var') || isempty(method)
    method = 'logmexp';
end
if any(size(Y)==1)
    Y = reshape(Y, [], 1);
else
    Y = Y';
end

%% estimate the noise
[psdx, ff] = pwelch(double(Y), [],[],[], 1);
indf = and(ff>=range_ff(1), ff<=range_ff(2));
switch method
    case 'mean'
        sn=sqrt(mean(psdx(indf, :)/2));
    case 'median'
        sn=sqrt(median(psdx(indf,:)/2));
    case 'logmexp'
        sn = sqrt(exp(mean(log(psdx(indf,:)/2))));    
    otherwise
        fprintf('wrong method! use logmexp instead.\n'); 
        sn = sqrt(exp(mean(log(psdx(indf,:)/2))));
end
sn = sn';
end

function [c, s, active_set] = oasisAR1(y, g, lam, smin, active_set)
%% Infer the most likely discretized spike train underlying an AR(1) fluorescence trace
% Solves the sparse non-negative deconvolution problem
%  min 1/2|c-y|^2 + lam |s|_1 subject to s_t = c_t-g c_{t-1} >=s_min or =0

%% inputs:
%   y:  T*1 vector, One dimensional array containing the fluorescence intensities
%withone entry per time-bin.
% OR %%
% len_active_set*4 matrix, active set

%   g:  scalar, Parameter of the AR(1) process that models the fluorescence ...
%impulse response.
%   lam:  scalar, sparsity penalty parameter lambda.
%   smin: scalar, optional, default 0
%miniumal non-zero activity within each bin (minimal 'spike size').
%   active_set: npool x 4 matrix, warm stared active sets

%% outputs
%   c: T*1 vector, the inferred denoised fluorescence signal at each time-bin.
%   s: T*1 vector, discetized deconvolved neural activity (spikes)
%   active_set: npool x 4 matrix, active sets

%% Authors: Pengcheng Zhou, Carnegie Mellon University, 2016
% ported from the Python implementation from Johannes Friedrich

%% References
% Friedrich J et.al., NIPS 2016, Fast Active Set Method for Online Spike Inference from Calcium Imaging

%% initialization
y = reshape(y, [], 1);
if isempty(y)
    T = sum(active_set(:,4)); 
else
T = length(y);
end
if ~exist('g', 'var') || isempty(g)
    g = estimate_time_constant(y);
elseif length(g)>1
    c = zeros(T,1); 
    s = zeros(T,1); 
    active_set = []; 
    return; 
end
if ~exist('lam', 'var') || isempty(lam);   lam = 0; end
if ~exist('smin', 'var') || isempty(smin);   smin = 0; end
if ~exist('active_set', 'var') || isempty(active_set)
    len_active_set = T;
    active_set = [y-lam*(1-g),ones(T,1),(1:T)',ones(T,1),(1:T)'-1, (1:T)'+1]; % each row is one pool: (vi, wi, t, l)
    active_set(end, :) = [y(end)-lam,1,T,1,T-1,nan] ;
    active_set(1,5) = nan;
else
    len_active_set = size(active_set,1);
    active_set(:,5) = [nan; (1:len_active_set-1)']; 
    active_set(:,6) = [(2:len_active_set)';nan]; 
end
idx = true(len_active_set,1);

%% run OASIS
ii = 1;
ii_next = active_set(ii,6);
while ~isnan(ii_next)
    % find the active set
    while (~isnan(ii_next)) && (active_set(ii_next,1)/active_set(ii_next,2)...
            >=active_set(ii,1)/active_set(ii,2)*g^(active_set(ii,4))+smin)
        active_set(ii_next,5) = ii;
        ii = ii_next; 
        ii_next = active_set(ii,6);
    end
    
    if isnan(ii_next); break; end
    
    %% merge pools
    active_set(ii,1) = active_set(ii,1) + active_set(ii_next,1)* (g^(active_set(ii,4)));
    active_set(ii,2) = active_set(ii,2) + active_set(ii_next,2)*(g^(2*active_set(ii,4)));
    active_set(ii,4) = active_set(ii,4) + active_set(ii_next,4);
    active_set(ii,6) = active_set(ii_next,6);
    idx(ii_next) = false;
    ii_next = active_set(ii,6);
    ii_prev = active_set(ii, 5);

    %% backtrack until violations fixed
    while (~isnan(ii_prev)) && (active_set(ii,1)/active_set(ii,2)<...
            max(0, active_set(ii_prev,1)/active_set(ii_prev,2)*g^(active_set(ii_prev,4)))+smin)
        ii_next = ii;
        ii = ii_prev;
        active_set(ii,1) = active_set(ii,1) + active_set(ii_next,1)* (g^(active_set(ii,4)));
        active_set(ii,2) = active_set(ii,2) + active_set(ii_next,2)*(g^(2*active_set(ii,4)));
        active_set(ii,4) = active_set(ii,4) + active_set(ii_next,4);
        active_set(ii,6) = active_set(ii_next,6);
        idx(ii_next) = false;
        
        ii_prev = active_set(ii, 5);
        ii_next = active_set(ii,6);
    end
end
active_set(~idx, :) = [];
len_active_set = size(active_set,1);

%% construct solution for all t
c = zeros(T, 1);
s = c;
for ii=1:len_active_set
    t0 = active_set(ii,3);
    tau = active_set(ii, 4);
    c(t0:(t0+tau-1)) = max(0,active_set(ii,1)/active_set(ii,2)) * (g.^(0:(tau-1)));
end

s(active_set(2:end,3)) = c(active_set(2:end,3)) - g*c(active_set(2:end,3)-1);
end


