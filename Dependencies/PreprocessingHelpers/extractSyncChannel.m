

function syncDat = extractSyncChannel(folder, numChans, syncChanIndex)
% extraChanIndices are 1-indexed

dataFiles = dir(fullfile(folder,'*.ap.*bin'));
if isempty(dataFiles)
    dataFiles = dir(fullfile(folder,'*.nidq.bin'));
end
for d = 1:length(dataFiles)
    filename = fullfile(dataFiles(d).folder,dataFiles(d).name);
    [folder,fn] = fileparts(filename);
    syncFname =  fullfile(folder, [fn '_sync.dat']);
    if exist(syncFname)
        fidOut = fopen(syncFname, 'r');
        fprintf(1,' loading %s\n', syncFname);
        syncDat = fread(fidOut, [1, Inf], 'int16=>int16'); % skipping other channels
        fclose(fidOut);
    else
        syncDat = extractSyncChannelFromFile(filename, numChans, syncChanIndex);
    end
    
end


disp(' done.')