function [fs,fc,prf,tau,bw,veff,ro,fdop] = ERSParameterExtractor(file)
% Open the parameter file to extract required parameters
fid = fopen(file,'r');

% Radar wavelength (satellite specific)
status = fseek(fid,720+500,'bof');
lambda = str2double(fread(fid,[1 16],'*char'));         % Wavelength (m)

% Pulse Repetition Frequency (satellite specific)
status = fseek(fid,720+934,'bof')|status;
prf = str2double(fread(fid,[1 16],'*char'));            % PRF (Hz)

% Range sampling rate (satellite specific)
status = fseek(fid,720+710,'bof')|status;
fs =str2double(fread(fid,[1 16],'*char'))*1e+06;        % Sampling Rate (Hz)

% Range Pulse length (satellite specific)
status = fseek(fid,720+742,'bof')|status;
tau = str2double(fread(fid,[1 16],'*char'))*1e-06;      % Pulse Width (sec)

% Range Gate Delay to first range cell
status = fseek(fid,720+1766,'bof')|status;
rangeGateDelay = str2double(fread(fid,[1 16],'*char'))*1e-03;   % Range Gate Delay (sec)

% Velocity X
status = fseek(fid,720+1886+452,'bof')|status;
xVelocity = str2double(fread(fid,[1 22],'*char'));    % xVelocity (m/sec)

% Velocity Y
status = fseek(fid,720+1886+474,'bof')|status;
yVelocity = str2double(fread(fid,[1 22],'*char'));    % yVelocity (m/sec)

% Velocity Z
status = fseek(fid,720+1886+496,'bof')|status;
zVelocity = str2double(fread(fid,[1 22],'*char'));    % zVelocity (m/sec)
fclose(fid);

% Checking for any file error
if(status==1)
    fs = NaN;
    fc = NaN;
    prf = NaN;
    tau = NaN;
    bw = NaN;
    veff = NaN;
    ro = NaN;
    fdop = NaN;
    return;
end

% Values specific to ERS satellites
slope = 4.19e+11;           % Slope of the transmitted chirp (Hz/s)
h = 790000;                 % Platform altitude above ground (m)
fdop = -1.349748e+02;       % Doppler frequency (Hz)

% Additional Parameters
Re = 6378144 ;              % Earth radius (m)

% Min distance
ro = time2range(rangeGateDelay);  % Min distance (m)  

% Effective velocity
v = sqrt(xVelocity^2+yVelocity^2+zVelocity^2);
veff = v*sqrt(Re/(Re+h));   % Effective velocity (m/sec)

% Chirp frequency
fc = wavelen2freq(lambda);  % Chirp frequency (Hz)     

% Chirp bandwidth
bw = slope*tau;             % Chirp bandwidth (Hz)
end