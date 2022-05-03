%% Automated
% Load all data
% Find available datasets (always using dates as folders)
clear DateOpt
for idx = 1:length(DataDir)
    DateOpt{idx} = cellfun(@(X) dir(fullfile(DataDir{idx},X,'*-*')),MiceOpt(DataDir2Use==idx),'UniformOutput',0);
end
DateOpt = cat(2,DateOpt{:});
DateOpt = cellfun(@(X) X([X.isdir]),DateOpt,'UniformOutput',0);
DateOpt = cellfun(@(X) {X.name},DateOpt,'UniformOutput',0);
RedoAfterClustering=1;
NewHistologyNeeded = 0; %Automatically to 1 after RedoAfterClustering

for midx = 1:length(MiceOpt)
    %% which probes?
    myKsDir = fullfile(LocalDir,MiceOpt{midx});
    subksdirs = dir(fullfile(myKsDir,'Probe*')); %This changed because now I suddenly had 2 probes per recording
    
    if strcmp(ProbeType{midx},'2_4S') && ~isempty(subksdirs) % These are my chronic mice, one dataset per mouse
        %% Loading data from kilosort/phy easily
        myKsDir = fullfile(LocalDir,MiceOpt{midx});
        subksdirs = dir(fullfile(myKsDir,'Probe*')); %This changed because now I suddenly had 2 probes per recording
        if length(subksdirs)<1
            clear subksdirs
            subksdirs.folder = myKsDir; %Should be a struct array
            subksdirs.name = 'Probe0';
        end
        for probeid = 1:length(subksdirs)
            myKsDir = fullfile(subksdirs(probeid).folder,subksdirs(probeid).name)
            if ~isdir(myKsDir)
                continue
            end
            
            %Saving directory
            thisprobe = subksdirs(probeid).name
            if (~NewHistologyNeeded && exist(fullfile(SaveDir,MiceOpt{midx},thisprobe,'HistoEphysAlignment.mat'))) && (~RedoAfterClustering || exist(fullfile(SaveDir,MiceOpt{midx},thisprobe,'CuratedResults.mat')))
                disp([MiceOpt{midx} ' already aligned... skip'])
                continue
            elseif RedoAfterClustering || NewHistologyNeeded
                myKsDir = fullfile(LocalDir,MiceOpt{midx},thisprobe);
                myClusFile = dir(fullfile(myKsDir,'cluster_info.tsv'));
                if isempty(myClusFile)
                    disp([MiceOpt{midx} ' is not yet curated with phy!!'])
                end
                NewHistologyNeeded = 1; %Automatically to 1 after RedoAfterClustering
            end
            
            %% Get cluster information
            clear params
            params.loadPCs=true;
            thisdate = [];
            PrepareClusInfo
            
            %% Get LFP?
            myLFDir = fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},'*','ephys');
            lfpD = dir(fullfile(myLFDir,'*','*','*.lf.*bin')); % ap file from spikeGLX specifically
            if isempty(lfpD)
                disp('No LFP data found, maybe it lives in .ap file?')
                lfpD = dir(fullfile(myLFDir,'*','*','*.ap.*bin')); % ap file from spikeGLX specifically
            end
            if isempty(lfpD)
                disp('No, really no data found..')
            elseif length(lfpD)>length(subksdirs)
                disp('Just take data from the last recording')
                lfpD = lfpD(end);
            elseif length(lfpD)<length(subksdirs)
                disp('Should be a different amount of probes?')
                keyboard
            else
                lfpD = lfpD(probeid);
            end
            
            %% Get Histology output
            GetHistologyOutput      %I created an extra function to have one line of code in the different scripts
            if ~histoflag
                disp([MiceOpt{midx} ' ' thisprobe 'No histology data, skip...'])
                continue
            end
        end
    else
        % For every date a different dataset
        Dates4Mouse = DateOpt{midx};
        for didx = 1:length(Dates4Mouse)
            % Within folders, look for 'RF mapping sessions'
            thisdate = Dates4Mouse{didx};
            %% Loading data from kilosort/phy easily
            myKsDir = fullfile(LocalDir,MiceOpt{midx},thisdate);
            subksdirs = dir(fullfile(myKsDir,'Probe*')); %This changed because now I suddenly had 2 probes per recording
            if length(subksdirs)<1
                clear subksdirs
                subksdirs.folder = myKsDir; %Should be a struct array
                subksdirs.name = 'Probe0';
            end
            for probeid = 1:length(subksdirs)
                myKsDir = fullfile(subksdirs(probeid).folder,subksdirs(probeid).name)
                if ~isdir(myKsDir)
                    continue
                end
                
                %Saving directory
                thisprobe = subksdirs(probeid).name
                if (~NewHistologyNeeded && exist(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,'HistoEphysAlignment.mat'))) && (~RedoAfterClustering || exist(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,'CuratedResults.mat')))
                    disp([MiceOpt{midx} ' ' thisdate ' already aligned... skip'])
                    continue
                elseif RedoAfterClustering || NewHistologyNeeded
                    myKsDir = fullfile(LocalDir,MiceOpt{midx},thisdate,thisprobe);
                    myClusFile = dir(fullfile(myKsDir,'cluster_info.tsv'));
                    if isempty(myClusFile)
                        disp([MiceOpt{midx} ' ' thisdate 'is not yet curated with phy!!'])
                    end
                    NewHistologyNeeded = 1; %Automatically to 1 after RedoAfterClustering
                end
                
                
                %% Get cluster information
                clear params
                params.loadPCs=true;
                params.thisdate = thisdate;
                PrepareClusInfo
                
                %% Get LFP?
                myLFDir = fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},thisdate,'ephys');
                lfpD = dir(fullfile(myLFDir,'*','*','*.lf.*bin')); % ap file from spikeGLX specifically
                if isempty(lfpD)
                    disp('No LFP data found')
                elseif length(lfpD)~=length(subksdirs)
                    disp('Should be a different amount of probes?')
                    keyboard
                else
                    lfpD = lfpD(probeid);
                end
                 
                %% Get Histology output
                GetHistologyOutput      %I created an extra function to have one line of code in the different scripts
                if ~histoflag
                    disp([MiceOpt{midx} ' ' thisdate ' ' thisprobe 'No histology data, skip...'])
                    continue
                end
            end
        end
    end
end
