function rawData = ERSDataExtractor(datafile,fs,fc,prf,tau,v,ro,doppler)
c = physconst('LightSpeed');                    % Speed of light

% Values specific to data file
totlines = 28652;                                % Total number of lines
numLines = 2048;                                 % Number of lines
numBytes = 11644;                                % Number of bytes of data
numHdr = 412;                                    % Header size
nValid = (numBytes-numHdr)/2 - round(tau*fs);    % Valid range samples

% Antenna length specific to ERS
L = 10;

% Calculate valid azimuth points
range = ro + (0:nValid-1) * (c/(2*fs));             % Computes range perpendicular to azimuth direction 
rdc = range/sqrt(1-(c*doppler/(fc*(2*v))^2));       % Squinted range 
azBeamwidth = rdc * (c/(fc*L)) * 0.8;               % Use only 80%  
azTau = azBeamwidth / v;                            % Azimuth pulse length 
nPtsAz = ceil(azTau(end) * prf);                    % Use the far range value
validAzPts = numLines - nPtsAz ;                    % Valid azimuth points  

% Start extracting
fid = fopen(datafile,'r');
status = fseek(fid,numBytes,'bof');     % Skipping first line
numPatch = floor(totlines/validAzPts);  % Number of patches           


if(status==-1)
   rawData = NaN;
   return;
end
rawData=zeros(numPatch*validAzPts,nValid);
% Patch data extraction starts
for patchi = 1:numPatch      
    fseek(fid,11644,'cof');
    data = fread(fid,[numBytes,numLines],'uint8')'; % Reading raw data file
    
    % Interpret as complex values and remove mean
    data = complex(data(:,numHdr+1:2:end),data(:,numHdr+2:2:end));
    data = data - mean(data(:));
    
    rawData((1:validAzPts)+((patchi-1)*validAzPts),:) = data(1:validAzPts,1:nValid);
    fseek(fid,numBytes + numBytes*validAzPts*patchi,'bof');
end
fclose(fid);
end