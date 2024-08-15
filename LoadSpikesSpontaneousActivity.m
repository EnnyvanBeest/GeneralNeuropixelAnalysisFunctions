clear DateOpt
DateOpt = arrayfun(@(X) dir(fullfile(DataDir{DataDir2Use(X)},MiceOpt{X},'*-*')),1:length(MiceOpt),'UniformOutput',0);

DateOpt = cellfun(@(X) X([X.isdir]),DateOpt,'UniformOutput',0);
DateOpt = cellfun(@(X) {X.name},DateOpt,'UniformOutput',0);
RedoAfterClustering=1;
Redo = 1; % Redo in general
NewHistologyNeeded = 0; %Automatically to 1 after RedoAfterClustering

%% Find RS sessions
%Predefine
pretrialtime = 2; %take up to x seconds prior trial
posttrialtime = 2; % take up to x seconds post trial
abortsession = 0;
for midx = 1:length(MiceOpt)
    Dates4Mouse = DateOpt{midx};
    for didx = 1:length(Dates4Mouse)
        % Within folders, look for 'RF mapping sessions'
        thisdate = Dates4Mouse{didx};
        subsess = dir(fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},Dates4Mouse{didx}));
        subsess(1:2) = []; %remove '.' and '..'
        flag = 0;
        ok = 0;
        RFsess = [];
        for sesidx=1:length(subsess)
            if strcmp(subsess(sesidx).name,'leftovers')
                continue
            end
            listfiles = dir(fullfile(subsess(sesidx).folder,subsess(sesidx).name));
            listfiles(~cell2mat(cellfun(@isempty,cellfun(@(X) strfind(X,'.pickle'),{listfiles(:).name},'UniformOutput',0),'UniformOutput',0))) = [];
            if any(~cell2mat(cellfun(@isempty,cellfun(@(X) strfind(X,'.p'),{listfiles(:).name},'UniformOutput',0),'UniformOutput',0)))
                idx = find(~cell2mat(cellfun(@isempty,cellfun(@(X) strfind(X,'.p'),{listfiles(:).name},'UniformOutput',0),'UniformOutput',0)));
                % read which x-file was used
                if length(idx)>1
                    disp('Too many files?!')
                    keyboard
                end
                fileID = fopen(fullfile(listfiles(idx).folder,listfiles(idx).name));
                A = fscanf(fileID,'%c');
                fclose(fileID)
                
                if any(strfind(A,'SparseNoise'))
                    RFsess = [RFsess {subsess(sesidx).name}];
                    flag = 1;
                    continue
                end
            else
                continue
            end
        end
        
        if ~flag
            continue
        end
        
        RSsess = [];
        for sesidx=1:length(RFsess)
            thisses=RFsess{sesidx};
            protocoldfile = dir(fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},thisdate,thisses,'Protocol.mat'));
            if isempty(protocoldfile)
                disp('No RF mapping file found... skip')
                continue
            end
            protocol = load(fullfile(protocoldfile.folder,protocoldfile.name));
            protocol = protocol.Protocol;
            
            %get number of row for c from parnames
            
            parrow = find(strcmp(protocol.parnames, 'c'));
            
            %check for 0 in row in pars
            dim = size(protocol.pars);
            sessions = dim(2);
            
            for session = 1:sessions
                if protocol.pars(parrow,session)==0
                    RSsess = [RSsess {thisses}];
                    ok = 1;
                end
            end
        end
        
        if ~ok
            continue
        end
        
        for sesidx = 1:length(RSsess)
            abortsession = 0;
            thisses = RSsess{sesidx}

            %% Timeline
            % Timeline to go to:
            clear Actualtime;
            timelinefile = dir(fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},thisdate,thisses,'*Timeline.mat'));
            timelineabsentwarning =0;
            if isempty(timelinefile)
                timelineabsentwarning =0;
                disp([MiceOpt{midx} ' ' thisdate   ' ' thisses ' has no timeline ... skipping'])
                continue
            end
            loadTimeline
            
            %% Align to Trial Onset times
            % Load from protocol trial info
            protocoldfile = dir(fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},thisdate,thisses,'Protocol.mat'));
            protocol = load(fullfile(protocoldfile.folder,protocoldfile.name));
            protocol = protocol.Protocol;
            ntrials = protocol.nrepeats; %This was not anymore the correct protocol so you have to reload it
            
            if isempty(contains(protocol.pardefs, 'Stimulus duration (s *10)'))
                disp('trial duration not found')
                keyboard
            else
                TrialDuration = protocol.pars(contains(protocol.pardefs, 'Stimulus duration (s *10)')); %in s*10
            end
            
            TrialDurations = cell(1, ntrials);
            TrialDurations(:) = {TrialDuration};
            TrialDurations = cell2mat(TrialDurations);
            %
            [starttrialidx,endtrialidx,trialid] = FindTrialOnAndOffSetsPhotoDiode(Actualtime,Timeline(:,ismember(AllInputs,'syncEcho')),TrialDurations);
            
            %% Checks
            if any((Actualtime(endtrialidx)-Actualtime(starttrialidx)) - TrialDurations'>1/tmSR*10)
                disp(['flips slightly drifting... Average of ' num2str(nanmean((Actualtime(endtrialidx)-Actualtime(starttrialidx)) - TrialDurations')) 'sec'])
            end
            
            if ntrials ~= length(starttrialidx)
                warning('Can''t find enough trials')
                keyboard
            end
            
            %% Loading data from kilosort/phy easily            
            myKsDir = fullfile(LocalDir,MiceOpt{midx},thisdate);
            subksdirs = dir(fullfile(myKsDir,'Probe*')); %This changed because now I suddenly had 2 probes per recording
            if length(subksdirs)<1
                clear subksdirs
                subksdirs.folder = myKsDir; %Should be a struct array
                subksdirs.name = 'Probe1';
            end
            for probeid = 1:length(subksdirs)
                myKsDir = fullfile(subksdirs(probeid).folder,subksdirs(probeid).name);
                
                if ~isdir(myKsDir)
                    continue
                end
                
                %Saving directory
                thisprobe = ['probe' num2str(probeid)]
                if ~Redo && exist(fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,thisprobe,'SpikeData.mat')) % replace with output
                    if ~RedoAfterClustering
                        continue
                    elseif RedoAfterClustering
                        myKsDir = fullfile(LocalDir,MiceOpt{midx},thisdate);
                        myClusFile = dir(fullfile(myKsDir,'cluster_info.tsv'));
                        if isempty(myClusFile)
                            disp('This data is not yet curated with phy!!')
                            continue
                        end
                        NewHistologyNeeded = 1; %Automatically to 1 after RedoAfterClustering
                    end
                end
                
                if ~isdir(fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,thisprobe))
                    mkdir(fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,thisprobe))
                end
                
                myLFDir = fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},thisdate,'ephys');
                sp = loadKSdir(myKsDir);
                % sp.st are spike times in seconds
                % sp.clu are cluster identities
                % spikes from clusters labeled "noise" have already been omitted
               
                %% Plotting a driftmap                
                [spikeTimes, spikeAmps, spikeDepths, spikeSites] = ksDriftmap(myKsDir);
                figure; plotDriftmap(spikeTimes, spikeAmps, spikeDepths);
                saveas(gcf,fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,thisprobe,'DriftMap.fig'))
                saveas(gcf,fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,thisprobe,'DriftMap.bmp'))
                
                %% basic quantification of spiking plot                
                depthBins = 0:40:3840;
                ampBins = 0:30:min(max(spikeAmps),800);
                recordingDur = sp.st(end);
                
                [pdfs, cdfs] = computeWFampsOverDepth(spikeAmps, spikeDepths, ampBins, depthBins, recordingDur);
                plotWFampCDFs(pdfs, cdfs, ampBins, depthBins);
                
                lfpFs = 2500;  % neuropixels phase3a - also 3b
                nChansInFile = 385;  % neuropixels phase3a, from spikeGLX
                
                %% Computing some useful details about spikes/neurons (like depths)                
                lfpD = dir(fullfile(myLFDir,'*','*','*.ap.bin')); % ap file from spikeGLX specifically
                if length(lfpD)~=length(subksdirs)
                    disp('Should be a different amount of probes?')
                    keyboard
                end
                lfpD = lfpD(probeid);
                %             if length(lfpD)>1
                %                 lfpD = lfpD(find(~cellfun(@isempty,cellfun(@(X) strfind(X,'Concat'),{lfpD(:).folder},'UniformOutput',0))));
                %             end
                %             if length(lfpD)<1
                %                 keyboard
                %             end
                [spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
                    templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);
                spikeCluster = sp.clu;
                
                %% Get cluster information                
                PrepareClusInfo

                %% Get Histology output
                GetHistologyOutput      %I created an extra function to have one line of code in the different scripts           
                if ~histoflag
                    disp([MiceOpt{midx} thisdate thisses thisprobe 'No histology data, skip...'])
                    continue
                end
                %% load synchronization data              
                SyncKSDataToTimeline % I created an extra function to have one line of code in the different scripts
                
                %if necessary
                if syncchanmissing % There were some sessions where the clock was outputted from imec and this signal was also written on flipper. Try to extract that, combined with neuronal data
                    syncchanmissingTrySyncAgain
                end
                if abortthissession
                    continue
                end
                if ~syncchanmissing && ~any(~isnan(spikeTimesCorrected))
                    warning('No Spikes in this session... continue')
                    continue
                end
              
                %remove NaNs
                spikeAmps(isnan(spikeTimesCorrected))=[];
                spikeCluster(isnan(spikeTimesCorrected)) = [];
                spikeDepths(isnan(spikeTimesCorrected))=[];
                spikeTimesCorrected(isnan(spikeTimesCorrected))=[];               
                
                %% Save Preprocessed data
                Trialtimes.starttime = Actualtime(starttrialidx); %Saving out trial start and endtimes helps you later (e.g. to define at what timepoint to start looking at the data)
                Trialtimes.endtime = Actualtime(endtrialidx);
                save(fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,thisprobe,'SpikeData.mat'),'tmSR','spikeTimesCorrected','Actualtime','Trialtimes','spikeCluster','Depth2AreaPerChannel','Depth2AreaPerUnit','clusinfo','spikeAmps','spikeDepths','-v7.3')
                
                
                %% Clear space
                clear sp
                clear SpikeRatePerPos
                clear SpikeRatePerTP
                clear spikeTimesCorrected
                clear allP
                clear ampBins
                clear cdfs
                clear clusinfo
                clear dataArray
                clear FlipperGLX
                clear Flippertimeline
                clear removevec
                clear rewardevents
                clear spikeAmps
                clear spikeCluster
                clear spikeDepths
                clear spikeSites
                clear spikesThisTrial
                clear spikeTimes
                clear spikeTimestmp
                clear SpikeTrialID
                clear SpikeTrialTime
                clear tempsUnW
                clear Actualtime
                clear allPowerEst
                clear allPowerVar
                clear F
                close all
                clear starttime
                clear h1
            end
        end
    end
    
end