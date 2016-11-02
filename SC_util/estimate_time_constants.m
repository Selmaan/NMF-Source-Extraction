function g = estimate_time_constants(y,p,sn,lags,nBinsAR,fR,minTau,maxTau)
    % estimate time constants from the sample autocovariance function
    
%% Default Parameters
    
if ~exist('nBinsAR','var') || isempty(nBinsAR)
    nBinsAR = 10; % # of bins for computing AR model coefficients
end

if ~exist('fR','var') || isempty(fR)
    fR = 6; % frame rate for conversion from s to frame tau
end

if ~exist('minTau','var') || isempty(minTau)
    minTau = 1; % in seconds
end

if ~exist('maxTau','var') || isempty(maxTau)
    maxTau = 4; % in seconds
end

minG = 1-1/(minTau*fR);
maxG = 1-1/(maxTau*fR);
%%
nI = floor(length(y)/nBinsAR);
yMat = reshape(y(1:nI*nBinsAR),nI,nBinsAR);
lags = lags + p;

for i=1:nBinsAR
    xc(i,:) = xcov(yMat(:,i),lags,'biased');
end
xc = median(xc)';

A = toeplitz(xc(lags+(1:lags)),xc(lags+(1:p))) - sn^2*eye(lags,p);
g = pinv(A)*xc(lags+2:end); 

if length(g) == 1
    g = max([g, minG]);
    g = min([g, maxG]);
end
