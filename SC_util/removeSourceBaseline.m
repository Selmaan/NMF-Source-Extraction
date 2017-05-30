function [fSub, f_] = removeSourceBaseline(f,binSize,basePrct)

nF = numel(f);
fOrig = f;
f(nF+1-rem(nF,binSize):end)=[];
f = reshape(f, binSize, []);
baseEstimate = prctile(f, basePrct, 1);

% 3rd order polynomial fits drift very well and is very robust (in contrast
% to exp()-based functions).
x = [(1:nF)', (1:nF)'.^2, (1:nF)'.^3];

xBin = x((binSize:binSize:length(baseEstimate)*binSize) - binSize/2, :);
warnState = warning('off', 'MATLAB:rankDeficientMatrix');

try
    b = robustfit(xBin, baseEstimate', 'bisquare', 2);
catch err
    keyboard,
end
warning(warnState);

f = fOrig;
f_ = [ones(nF, 1), x] * b;
f_ = f_';
fSub = f - f_;
end