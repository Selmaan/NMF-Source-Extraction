function [P,Y] = preprocess_data(Y,p,options)

% data pre-processing for:
% (ii)  identifying saturated pixels
% (iii) estimating noise level for every pixel

% Author: Eftychios A. Pnevmatikakis
%           Simons Foundation, 2015

defoptions.noise_range = [0.33,0.5];            % frequency range over which to estimate the noise
defoptions.noise_method = 'logmexp';            % method for which to estimate the noise level
defoptions.block_size = [64,64];
defoptions.split_data = 0;

if nargin < 3 || isempty(options); options = defoptions; end
if nargin < 2 || isempty(p); p = 0; end
P.p = p;

if ~isfield(options,'noise_range'); options.noise_range = defoptions.noise_range; end
if ~isfield(options,'noise_method'); options.noise_method = defoptions.noise_method; end
if ~isfield(options,'block_size'); options.block_size = defoptions.block_size; end
if ~isfield(options,'split_data'); split_data = defoptions.split_data; else split_data = options.split_data; end

%% indentify saturated pixels

P.pixels = find_unsaturatedPixels(Y);                % pixels that do not exhibit saturation

%% estimate noise levels
[sn,psx] = get_noise_fft(Y,options);
P.sn = sn(:);