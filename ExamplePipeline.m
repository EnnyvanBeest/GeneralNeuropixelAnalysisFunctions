
%% Predefine
DataDir = [path/to/data];
SaveDir = [path/to/outputfolder];
pretrialtime = 2; %take up to x seconds prior trial
posttrialtime = 2; % take up to x seconds post trial

%% Loading data from kilosort/phy easily
myKsDir = fullfile(DataDir,SessionName);

%Saving directory
if ~isdir(fullfile(SaveDir)
    mkdir(fullfile(SaveDir)
end

myLFDir = fullfile(DataDir,'ephys');
sp = loadKSdir(myKsDir); %Function from spikes toolbox
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
lfpD = dir(fullfile(DataDir,'*','*','*.ap.bin')); % ap file from spikeGLX specifically

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
    disp('No histology data, skip...')
    continue
end
%% load synchronization data
timelinefile = dir(fullfile(DataDir,'*Timeline.mat')); %Some file with all timestamps of different signals. E.g. actual time, photodiode input, encoder input, etc.
SyncKSDataToTimeline 

if ~syncchanmissing && ~any(~isnan(spikeTimesCorrected))
    warning('No Spikes in this session... continue')
    continue
end

%% Align to Trial Onset times
% Load from protocol trial info
protocoldfile = dir(DataDir,'Protocol.mat'); %Some experimental protocol file
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
[starttrialidx,endtrialidx,trialid] = FindTrialOnAndOffSetsPhotoDiode(Actualtime,Timeline(:,ismember(AllInputs,'photodiode')),TrialDurations);

%%
if any((Actualtime(endtrialidx)-Actualtime(starttrialidx)) - TrialDurations'>1/tmSR*10)
    disp(['flips slightly drifting... Average of ' num2str(nanmean((Actualtime(endtrialidx)-Actualtime(starttrialidx)) - TrialDurations')) 'sec'])
end

if ntrials ~= length(starttrialidx)
    error('Can''t find enough trials')
end

if syncchanmissing % There were some sessions where the clock was outputted from imec and this signal was also written on flipper. Try to extract that, combined with neuronal data
    
    % Aproximate start of session relative to ephys
    %             fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,'SpikeData.mat')
    allsess = dir(fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},thisdate)); %If you use continuous ephys but have separate experimental sessions 
    sesnrs = (cellfun(@(X) str2num(X),{allsess(:).name},'UniformOutput',0));
    sesnrs(cellfun(@isempty,sesnrs))={[1000]};
    [~,idx] = sort([sesnrs{:}]);
    allsess = allsess(idx);
    totaldur = 0;
    for ij = 1:length(allsess)
        tmptimeline = dir(fullfile(allsess(ij).folder,allsess(ij).name,'*_Timeline.mat'));
        if ~isempty(tmptimeline)
            tmptl = load(fullfile(tmptimeline(1).folder,tmptimeline(1).name));
            if strcmp(allsess(ij).name,thisses)
                break
            end
            totaldur = totaldur+tmptl.Timeline.lastTimestamp;
        end
    end
    
    % Get clock signal that leaked into flipper
    Flippertimelineunfiltered = Timeline(:,ismember(AllInputs,'flipper'));
    imecclock = Flippertimelineunfiltered-Flippertimeline.*5;
    tmpmean = nanmedian(imecclock(:));
    if tmpmean<0
        imecclock(imecclock>tmpmean)=1;
        imecclock(imecclock<tmpmean)=0;
    else
        imecclock(imecclock<tmpmean)=0;
        imecclock(imecclock>tmpmean)=1;
    end
    
    imecclock = lowpass(imecclock,2,tmSR);
    tmpmean = nanmedian(imecclock(:));
    imecclock(imecclock<tmpmean)=0;
    imecclock(imecclock>tmpmean)=1;
    % Bin spiketimes in actual time
    timeBinSize = 1/tmSR; % seconds
    
    % GLM to predict spike times (timeXUnit) based on timeline
    % events --> to find best matching signal to timeline in
    % absence of flipper signal
    itioff = zeros(1,length(Actualtime));
    for trid=1:ntrials
        itioff(starttrialidx(trid):endtrialidx(trid))=1;
    end
    
    %Rewards
    RewardSignal = (Timeline(:,ismember(AllInputs,'rewardEcho')));
    RewardSignal(RewardSignal<3)=0;
    RewardSignal(RewardSignal>=3)=1;
    
    %RunningSpeed
    RotarySignal = (Timeline(:,ismember(AllInputs,'rotaryEncoder')));
    Speed = nan(1,length(RotarySignal));
    windowsz = tmSR; %per second
    parfor tp = 1:length(RotarySignal)
        tptotake = [tp-windowsz./2:tp+windowsz./2];
        tptotake(tptotake<1|tptotake>length(RotarySignal))=[];
        Speed(tp) = nanmean(diff(RotarySignal(uint16(tptotake))));
    end
    Speed = (Speed-nanmin(Speed(:)))./(nanmax(Speed(:))-nanmin(Speed(:))); %Normalize
    
    
    % all spikes allegedly in this session
    extratime = 0;
    tmpspikeshere= spikeTimestmp(spikeTimestmp>totaldur-extratime&spikeTimestmp<totaldur-extratime+tmptl.Timeline.lastTimestamp);
    clusterIDhere = spikeCluster(spikeTimestmp>totaldur-extratime&spikeTimestmp<totaldur-extratime+tmptl.Timeline.lastTimestamp);
    %             visid = find(cell2mat(cellfun(@isempty,cellfun(@(X) strfind(X,'vis'),Depth2AreaPerUnit.Area,'UniformOutput',0),'UniformOutput',0)));
    %             visid(~ismember(visid,Good_ID))=[];
    %             tmpspikeshere = tmpspikeshere(ismember(clusterIDhere,visid));
    %             clusterIDhere = clusterIDhere(ismember(clusterIDhere,visid));
    
    %spikecount per unit
    unitid = unique(clusterIDhere);
    spikecount = zeros(length([min(tmpspikeshere)-timeBinSize:timeBinSize:max(tmpspikeshere)+timeBinSize]),length(unitid));
    for i = 1:length(unitid)
        tmp = [0 (histcounts(tmpspikeshere(clusterIDhere==unitid(i)),[min(tmpspikeshere)-timeBinSize:timeBinSize:max(tmpspikeshere)+timeBinSize]))];
        spikecount(tmp>0,i) = 1;
    end
    
    
    %find maximum correlation between rotary signal and average spikecount
    [C21,lag21] = xcorr(Speed,nanmean(spikecount,2));
    C21 = C21./nanmax(C21);
    [maxval,maxid] = nanmax(C21);
    lag21(maxid)
    
    %find maximum correlation between itioff signal and average spikecount
    [C31,lag31] = xcorr(itioff,nanmean(spikecount,2));
    C31 = C31./nanmax(C31);
    [maxval,maxid] = nanmax(C31);
    lag31(maxid)
    
    %find maximum correlation between itioff signal and average spikecount
    [C41,lag41] = xcorr(RewardSignal',nanmean(spikecount,2));
    C41 = C41./nanmax(C41);
    [maxval,maxid] = nanmax(C41);
    lag41(maxid)
    
    % Sum together
    Cortot = C21+C31+C41;
    [maxval,maxid] = nanmax(Cortot);
    extratime = lag41(maxid)./tmSR; %in seconds
    
    figure;
    plot(lag41,Cortot)
    
    figure;
    disp(['Found session to start ' num2str(totaldur+extratime) ' seconds into recording'])
    delayestimate = round((totaldur+extratime)*tmSR);
    
    for i = 1:per2check:max(ceil(Actualtime))
        indx=find(Actualtime>=i-1 & Actualtime<i+per2check-1);
        indx(indx>length(Actualtime))=[];
        indx(indx>length(syncDatImec))=[];
        
        if isempty(indx)
            break
        end
        
        if exist('h1')
            delete(h1)
            delete(h2)
        end
        % Check every x timepoints what the difference is between the
        % two flippers
        indx2 = indx+delayestimate;
        indx(indx2<1)=[];
        indx2(indx2<1)=[];
        indx(indx2>length(syncDatImec))=[];
        indx2(indx2>length(syncDatImec))=[];
        if isempty(indx2)
            continue
        end
        
        delayindx = finddelay(imecclock(indx),syncDatImec(indx2));
        
        %D = FINDDELAY(X,Y), where X and Y are row or column vectors of length
        %   LX and LY, respectively, returns an estimate of the delay D between X
        %   and Y, where X serves as the reference vector. If Y is delayed with
        %   respect to X, D is positive. If Y is advanced with respect to X, D is
        %   negative.
        
        indx2 = indx2+delayindx;
        indx(indx2<1)=[];
        indx2(indx2<1)=[];
        indx(indx2>length(syncDatImec))=[];
        indx2(indx2>length(syncDatImec))=[];
        if isempty(indx2)
            continue
        end
        
        h1 = plot(Actualtime(indx),imecclock(indx),'b-');
        hold on
        h2 = plot(Actualtime(indx),syncDatImec(indx2)+2,'r-');
        %sanity check
        if warningflag
            AcceptableYN = '';
            while ~ismember(AcceptableYN,{'y','n'})
                AcceptableYN = input('Is this alignment acceptable? (y/n)','s');
                if strcmpi(AcceptableYN,'n')
                    disp('Not acceptable, skipping session')
                    keyboard
                    abortthissession = 1;
                    break
                elseif strcmp(AcceptableYN,'y')
                    disp('Acceptable, continuing...')
                    warningflag = 0;
                end
            end
        end
        if abortthissession
            break
        end
        
        %Convert back to time in imec data
        TL2ImecTime = (delayindx+delayestimate)./tmSR;
        Imec2TLTime = (-delayindx-delayestimate)./tmSR;
        
        %Find spikes in IMEC time, that fall in this window in
        %TIMELINE space
        spikeindx = find(spikeTimestmp>=TL2ImecTime+(i-1)&spikeTimestmp<=TL2ImecTime+(i+per2check-1));
        if ~isempty(spikeindx)
            spikeTimesCorrected(spikeindx)=spikeTimestmp(spikeindx)+Imec2TLTime;
        end
    end
end

if abortthissession
    continue
end

nclus = length(Good_ID);

%remove NaNs
spikeAmps(isnan(spikeTimesCorrected))=[];
spikeCluster(isnan(spikeTimesCorrected)) = [];
spikeDepths(isnan(spikeTimesCorrected))=[];
spikeTimesCorrected(isnan(spikeTimesCorrected))=[];

%% Save Preprocessed data
Trialtimes.starttime = Actualtime(starttrialidx); %Saving out trial start and endtimes helps you later (e.g. to define at what timepoint to start looking at the data)
Trialtimes.endtime = Actualtime(endtrialidx);
save(fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,thisprobe,'SpikeData.mat'),'tmSR','spikeTimesCorrected','Actualtime','Trialtimes','spikeCluster','Depth2AreaPerChannel','Depth2AreaPerUnit','clusinfo','-v7.3')

