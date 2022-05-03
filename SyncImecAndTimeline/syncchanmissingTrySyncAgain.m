% syncchanmissingTrySyncAgain



% Aproximate start of session relative to ephys
%             fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,'SpikeData.mat')
allsess = dir(fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},thisdate,thisprobe));
if isempty(allsess)
    allsess = dir(fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},thisdate));
end
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
RotarySignalPerSec = nan(1,length(RotarySignal));
windowsz = tmSR; %per second
parfor tp = 1:length(RotarySignal)
    tptotake = [tp-windowsz./2:tp+windowsz./2];
    tptotake(tptotake<1|tptotake>length(RotarySignal))=[];
    RotarySignalPerSec(tp) = nanmean(diff(RotarySignal(uint16(tptotake))));
end
RotarySignalPerSec = (RotarySignalPerSec-nanmin(RotarySignalPerSec(:)))./(nanmax(RotarySignalPerSec(:))-nanmin(RotarySignalPerSec(:)));

% all spikes allegedly in this session
extratime = 0;
tmpspikeshere= spikeTimestmp(spikeTimestmp>totaldur-extratime&spikeTimestmp<totaldur-extratime+tmptl.Timeline.lastTimestamp);
clusterIDhere = spikeCluster(spikeTimestmp>totaldur-extratime&spikeTimestmp<totaldur-extratime+tmptl.Timeline.lastTimestamp);
%             visid = find(cell2mat(cellfun(@isempty,cellfun(@(X) strfind(X,'vis'),Depth2AreaPerUnit.Area,'UniformOutput',0),'UniformOutput',0)));
%             visid(~ismember(visid,Good_ID))=[];
%             tmpspikeshere = tmpspikeshere(ismember(clusterIDhere,visid));
%             clusterIDhere = clusterIDhere(ismember(clusterIDhere,visid));

%spikecount per unit
spikecount = zeros(length([min(tmpspikeshere)-timeBinSize:timeBinSize:max(tmpspikeshere)+timeBinSize]),length(Good_IDx));
for i = 1:length(Good_IDx)
    tmp = [0 (histcounts(tmpspikeshere(clusterIDhere==Good_IDx(i)),[min(tmpspikeshere)-timeBinSize:timeBinSize:max(tmpspikeshere)+timeBinSize]))];
    spikecount(tmp>0,i) = 1;
end

figure;

%find maximum correlation between rotary signal and average spikecount
if any(RotarySignalPerSec)
    [C21,lag21] = xcorr(RotarySignalPerSec,nanmean(spikecount,2));
    C21 = C21./nanmax(C21);
    [maxval,maxid] = nanmax(C21);
    lag21(maxid)
    subplot(4,1,1)
    plot(lag21,C21)
    title('Rotary vs. Spikes')
else
    C21=nan(1,length(spikecount)*2-1);
end

%find maximum correlation between itioff signal and average spikecount
if any(itioff)
    [C31,lag31] = xcorr(itioff,nanmean(spikecount,2));
    C31 = C31./nanmax(C31);
    [maxval,maxid] = nanmax(C31);
    lag31(maxid)
    subplot(4,1,2)
    plot(lag31,C31)
    title('ITItime vs. Spikes')
else
    C31=nan(1,length(spikecount)*2-1);

end

%find maximum correlation between itioff signal and average spikecount
if any(RewardSignal)
    [C41,lag41] = xcorr(RewardSignal',nanmean(spikecount,2));
    C41 = C41./nanmax(C41);
    [maxval,maxid] = nanmax(C41);
    lag41(maxid)
    
    subplot(4,1,3)
    plot(lag41,C41)
    title('RewardSignal vs. Spikes')
else
    C41=nan(1,length(C21));
end

% Sum together
Cortot = nansum(cat(1,C21,C31,C41),1);
[maxval,maxid] = nanmax(Cortot);
extratime = lag31(maxid)./tmSR; %in seconds
subplot(4,1,4)
plot(lag31,Cortot)
title('Total Correlation SpikeData with Timeline Data versus lag')

for subid=1:4
    subplot(4,1,subid)
    hold on
    line([extratime extratime],get(gca,'ylim'))
end

figure;
disp(['Found session to start ' num2str(totaldur-extratime) ' seconds into recording'])
delayestimate = round((totaldur-extratime)*tmSR);

%% Correlation on whole set before fine tuning
indx=1:length(Actualtime);
indx(indx>length(Actualtime))=[];
indx(indx>length(syncDatImec))=[];

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


delayindx = finddelay(imecclock(indx),syncDatImec(indx2));
delayestimate = delayindx+delayestimate;
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

if corr(double(imecclock(single(indx))),double(syncDatImec(single(indx2))'))>0.75
    warningflag=0
    disp('Correlation okay.. continue')
elseif corr(double(imecclock(single(indx))),double(syncDatImec(single(indx2))'))<0.25
    disp('Correlation bad, no need to check manually, skip..')
    abortthissession = 1;
    return
else
    disp('Correlation not great...')
end

%% Fine tuning
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
    if warningflag && i>10
        drawnow
        AcceptableYN = '';
        while ~ismember(AcceptableYN,{'y','n'})
            AcceptableYN = input('Is this alignment acceptable? (y/n)','s');
            if strcmpi(AcceptableYN,'n')
                disp('Not acceptable, skipping session')
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
