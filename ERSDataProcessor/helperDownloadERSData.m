function helperDownloadERSData(outputFolder,DataURL)
% Download the data set from the given URL to the output folder.

    radarDataZipFile = fullfile(outputFolder,'ERSData.zip');
    
    if ~exist(radarDataZipFile,'file')
        
        disp('Downloading ERS data (134 MiB)...');
        websave(radarDataZipFile,DataURL);
        unzip(radarDataZipFile,outputFolder);
    end
end