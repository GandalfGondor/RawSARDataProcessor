function slcimg = myRangeMigrationLFM(raw,waveform,fc,v,rc,varargin)
%rangeMigrationLFM  Range migration image formation algorithm for LFM
%   SLCIMG = rangeMigrationLFM(RAW,WAVEFORM,FC,V,RC) generates a
%   single-look complex (SLC) image of raw synthetic aperture radar (SAR)
%   data using the range migration algorithm. The input RAW contains the
%   unfocused IQ raw data collected by a synthetic aperture system. RAW is
%   of size M-by-N. M is the number of along-range samples and N is the
%   number of pulses received as the platform moves along the cross-range
%   direction. M must be greater than 1. The input WAVEFORM is a
%   phased.LinearFMWaveform object. The input FC is a positive real scalar
%   that specifies the operating frequency in Hz. The input V is positive
%   real scalar that specifies the platform velocity in meters per second.
%   The input RC is a positive real scalar that specifies the distance in
%   meters between the beam center on the ground and the radar. The output
%   SLCIMG is a matrix of the same size as RAW that contains the focused
%   data processed by the range migration algorithm. The output SLCIMG
%   represents an SLC image.
%
%   SLCIMG = rangeMigrationLFM(...,'SquintAngle',SQA) specifies the squint
%   angle of the antenna from the broadside direction as a real scalar in
%   degrees. The default value of SQA is 0 degrees, which corresponds to a
%   broadside antenna. SQA takes values in the range (-90,90).
%
%   SLCIMG = rangeMigrationLFM(...,'PropagationSpeed',C) specifies the
%   signal propagation speed C in meters per second as a positive scalar.
%   The default value of C is the speed of light.
%
%   % Example:
%   %   Generate single-look complex image for the simulated unfocused SAR
%   %   raw data.
%
%   load('RangeMigrationLFMExampleData.mat')
%   slcimg = rangeMigrationLFM(raw,waveform,fc,v,rc);
%
%   %   Generate single-look complex image
%   imagesc(abs(slcimg))
%   title('SLC Image')
%   xlabel('Cross-Range Samples')
%   ylabel('Range Samples')
%
%   See also rangeMigrationSFM, rangeMigrationFMCW.

%   Copyright 2021 The Mathworks Inc.

%   References
%   [1] Cumming, Ian G., and Frank Hay-chee Wong. "Digital Processing of
%       Synthetic Aperture Radar Data: Algorithms and Implementation."
%       Artech House Remote Sensing Library. Boston: Artech House, 2005.
%   [2] Tolman, Matthew A. "A Detailed Look at the Omega-k Algorithm for
%       Processing Synthetic Aperture Radar Data." Master's Thesis, Brigham
%       Young University, 2008.

%#codegen

narginchk(5,9);

[sqa, c, expParam] = parseInputs(raw,waveform,fc,v,rc,varargin{:}); % Validate required inputs.

fs = waveform.SampleRate;           % Sampling frequency
numrange = size(raw,1);  % Number of along range samples
numcrossrng = size(raw,2);    % Number of cross-range samples
imagedata = zeros(numrange,numcrossrng,'like',1+1i);
frstart = fc-fs/2;
frend = fc+fs/2-fs/numrange;

prf = waveform.PRF;                 % Pulse repetition frequency

% Pulse compression
ref = calcMatchedFilterCoefficient(waveform);
for i = 1:numcrossrng
    convdata = conv(raw(:,i),ref);
    imagedata(:,i) = convdata(length(ref):length(convdata));
end

%% Squint angle compensation
kc = (2*pi*fc)/c;
slowtime = ((0:(numcrossrng-1))'/prf);
slowtime = slowtime-mean(slowtime);
imagedata = imagedata.*exp(-1i.*2*(kc)*sin(deg2rad(sqa))*repmat(v*slowtime,1,numrange)).';

rangeFocusProcessing(imagedata, expParam);
% save('rangeFocus.mat', 'imagedata');

% 替换 imagedata, 仅处理包含目标的。
imagedata = replaceSelectAreas(imagedata, expParam.selectArea);

figure;
imagesc(abs(imagedata));
axis image;colormap('gray');colorbar;
title('Range Compress Image');
ylabel('Range bin');
xlabel('Azimuth bin');
savefig(sprintf('./%s/result/focusTarget.fig', expParam.fileName));

% 
%% Converting to frequency domain
imagedata = fftshift(fft(imagedata,[],1),1);
imagedata = fftshift(fft(imagedata,[],2),2);

%% Range frequency span
frequencyRange = linspace(frstart,frend,numrange).';
krange = 2*(2*pi*frequencyRange)/c;

%% Wavenumber processing
kus = 2*(kc)*sin(deg2rad(sqa));
kcrossrng = 2*pi*linspace(-prf/2,prf/2-prf/numcrossrng,numcrossrng)./v;
krange = repmat(krange,1,numcrossrng);
kcrossrng = repmat(kcrossrng,numrange,1);
kx = krange.^2-(kcrossrng+kus).^2;
kx = sqrt(kx.*(kx > 0));
negativeRowIndices = any(kx <= 0);
kFinal = exp(1i*(kx.*cos(deg2rad(sqa))+(kcrossrng).*sin(deg2rad(sqa))).*rc);
kfin = kx.*cos(deg2rad(sqa))+(kcrossrng+kus).*sin(deg2rad(sqa));

%% Cross-range Compression
imagedata = (imagedata).*kFinal;

%% Interpolation
for i = 1:size((imagedata),2)
    if negativeRowIndices(i) == false
        imagedata(:,i) = interp1(kfin(:,i),imagedata(:,i),krange(:,1));
    end
end

imagedata(isnan(imagedata)) = 1e-30; % Replace NaN with small number
imagedata = imagedata.*exp(-1i*(krange).*rc);
slcimg = (ifft2(imagedata));
end

%--------------------------------------------------------------------------
function [sqa, c, expParam] = parseInputs(raw,waveform,fc,v,rc,varargin)

funName = 'rangeMigrationLFM';

validateattributes(raw, {'double','single'}, {'nonempty', '2d', 'nonnan',...
    'finite'},funName, 'RAW',1);

validateattributes(waveform,{'phased.LinearFMWaveform'},{},funName,...
    'WAVEFORM',2);

validateattributes(waveform.PRF,{'double'},{'scalar'});

validateattributes(fc,{'double'},{'nonempty','scalar','real',...
    'positive','nonnan','finite'},funName,'FC',3);

validateattributes(v,{'double'},{'nonempty','scalar','real',...
    'positive','nonnan','finite'},funName,'V',4);

validateattributes(rc,{'double'},{'nonempty','scalar','real',...
    'positive','nonnan','finite'},funName,'RC',5);

% Define default values for optional inputs.
default.SquintAngle = 0;
default.PropagationSpeed = physconst('LightSpeed');

parseoptions = struct('CaseSensitivity',false, ...
    'PartialMatching','none', ...
    'StructExpand',false, ...
    'IgnoreNulls',true,...
    'SupportOverrides',false);
opArgs = {};
params = {'SquintAngle','ExpParam','PropagationSpeed'};

pstruct = coder.internal.parseInputs(opArgs,params,parseoptions,varargin{:});
results = coder.internal.vararginToStruct(pstruct,default,varargin{:});

sqa = results.SquintAngle;
c = results.PropagationSpeed;
expParam = results.ExpParam;

validateattributes(sqa, {'double'}, {'nonempty', 'scalar', 'real', ...
    'nonnan','finite','>', -90, '<', 90}, funName, 'SQA');

validateattributes(c,{'double'},{'nonempty','scalar','real',...
    'positive','nonnan','finite'},funName,'C');

end