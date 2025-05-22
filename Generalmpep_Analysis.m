%% User Input
% Load all data
OriSetting = PipelineParams;

% Find available datasets (always using dates as folders)
RedoAfterClustering=0;
Redo = 1; % Redo in general
RedoTable = 1;
NewHistologyNeeded = 0; %Automatically to 1 after RedoAfterClustering
%Predefine
SaveRFDir = SaveDir
abortsession = 0;
timeBinSize = 1/100;
pretrialtime = 0.2; %take up to x seconds prior trial
posttrialtime = 0.2; % take up to x seconds post trial
% Build wavelet filter bank
freqlims = [0.85 9.79];
nv = 48;
plotthis = 0; %For intermediate step plots (only for debugging)
CopyToXFileSpecific = 0;
%% Automated
% Build filterbank for wavelet transforms
clear DateOpt
DateOpt = arrayfun(@(X) dir(fullfile(DataDir{DataDir2Use(X)},MiceOpt{X},'*-*')),1:length(MiceOpt),'UniformOutput',0);
DateOpt = cellfun(@(X) X([X.isdir]),DateOpt,'UniformOutput',0);
DateOpt = cellfun(@(X) {X.name},DateOpt,'UniformOutput',0);
for midx = 1:length(MiceOpt)
    Dates4Mouse = DateOpt{midx};
    for didx = 1:length(Dates4Mouse)
        % Within folders, look for 'RF mapping sessions'
        thisdate = Dates4Mouse{didx}
        subsess = dir(fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},Dates4Mouse{didx}));
        subsess(1:2) = []; %remove '.' and '..'
        flag = 0;
        mpepsess = [];
        for sesidx=1:length(subsess)
            if strcmp(subsess(sesidx).name,'leftovers')
                continue
            end
            listfiles = dir(fullfile(subsess(sesidx).folder,subsess(sesidx).name,'*.p')); %Find protocol file
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

                if ~any(strfind(A,'SparseNoise'))
                    mpepsess = [mpepsess {subsess(sesidx).name}];
                    flag = 1;
                    continue
                else
                    disp('RF session, use different script... continue')
                end
            else
                continue
            end
        end

        if ~flag
            continue
        end

        for sesidx = 1:length(mpepsess)
            close all
            abortsession = 0;
            thisses = mpepsess{sesidx}

            %% Timeline
            % Timeline to go to:
            timelinefile = dir(fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},thisdate,thisses,'*Timeline.mat'));
            timelineabsentwarning =0;
            if isempty(timelinefile)
                timelineabsentwarning =0;
                disp([MiceOpt{midx} ' ' thisdate   ' ' thisses ' has no timeline ... skipping'])
                continue
            end
            loadTimeline

            % Microphone data
            [MicData,MicFS,MicnBits,SP] = LoadMicData(fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},thisdate,thisses),0,0);
            if ~isempty(SP)
                f=SP.F;
                Pw=abs(SP.S);
                t = abs(SP.T);
            else
                f=[];
                Pw = [];
                t = [];
            end

            %% Load Protocol and identify unique conditions
            Protocol = load(fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},thisdate,thisses,'Protocol.mat'));
            Protocol = Protocol.Protocol;

            % extract trial information
            ntrials = prod(size(Protocol.seqnums));
            ParFullnames = Protocol.pardefs;
            ParDefs = Protocol.parnames;
            Pars = Protocol.pars;
            ncond = size(Pars,2);
            xfile = Protocol.xfile;

            %% For sound, check whether these frequencies can be found
            if any(ismember(ParDefs,'SoundFreq')) && ~isempty(Pw)
                FreqsUsed = unique(Pars(ismember(ParDefs,'SoundFreq'),:))./10*1000; %in Hz
                avgpw = nanmean(Pw,2);
                avgdf = [0; diff(avgpw)];
                figure; plot(f,nanmean(Pw,2));
                hold on;
                peaks = [];
                id=find(avgdf>nanmean(avgdf)+1.5*nanstd(avgdf),1,'first');
                while id <length(avgdf)
                    endidx = find(avgdf(id:end)<0,1,'first')+id;
                    [~,pkidx] = max(avgdf(id:endidx));
                    peaks = [peaks f(pkidx+id)];
                    id = endidx+find(avgdf(endidx:end)>nanmean(avgdf(endidx:end))+1.5*nanstd(avgdf(endidx:end)),1,'first');;
                end
                plot(peaks,avgpw(ismember(f,peaks)),'r*')
                [r,c] = find(abs(FreqsUsed-peaks')<1000);
                peaks = peaks(unique(r));
                plot(peaks,avgpw(ismember(f,peaks)),'g*')
                title('Microphone data')
                xlabel('frequency (Hz)')
                ylabel('|P|')
                makepretty

                if isempty(peaks)
                    disp('There was no sound!!')
                    Pars(ismember(ParDefs,'SoundFreq'),:)=0;
                    Pars(ismember(ParDefs,'ModFreq'),:)=0;
                end
            end

            %%
            sequence = Protocol.seqnums;
            [condindx,repidx] = arrayfun(@(X) find(sequence==X),1:ntrials,'UniformOutput',0);
            nrep = max(cell2mat(repidx));
            % Factors of relevance
            parrel = [];
            RelName = {};
            for parid = 1:size(Pars,1)
                if length(unique(Pars(parid,:)))>1
                    parrel = [parrel parid];
                    RelName = {RelName{:} ParDefs{parid}};
                end
            end
            notrelid = [];
            [CondOpt,idx,cdx] = arrayfun(@(X) unique(Pars(X,:),'stable'),parrel,'UniformOutput',0);
            for parid=1:length(parrel)
                for parid2=1:length(parrel)

                    if ~any((cdx{parid} == cdx{parid2})==0)
                        if nanvar(CondOpt{parid})<nanvar(CondOpt{parid2}) %Categorical, remove this one
                            notrelid=[notrelid parrel(parid)]; %Remove duplicate conditions (some conditions are 1:1 match, like amplitude and frequency for sound)
                        end
                    end
                end
            end
            if any(notrelid)
                RelName(ismember(parrel,notrelid))=[];
                parrel(ismember(parrel,notrelid))=[];

            end

            %Remove dcycl - not ideal but for now
            parrel(ismember(RelName,'dcyc'))=[];
            RelName(ismember(RelName,'dcyc'))=[];

            condindx = cell2mat(condindx); %Condition index in trial ordedr
            TrialDurations = Pars(strcmp(ParDefs,'dur'),condindx)/10;

            ControlCond = find(nanmean(Pars(ismember(ParDefs,{'cr','cb','cg'}),:),1)==0 | nanmean(Pars(ismember(ParDefs,{'lr','lb','lg'}),:),1)==0 ); %^Contrast or luminance 0
            RelName(ismember(parrel,find(ismember(ParDefs,{'cr','cb','cg','lr','lb','lg'}))))=[];
            parrel(ismember(parrel,find(ismember(ParDefs,{'cr','cb','cg','lr','lb','lg'}))))=[];
            newtimevec = -pretrialtime:timeBinSize:max(TrialDurations)+posttrialtime;
            timeEdges = -pretrialtime-timeBinSize/2:timeBinSize:max(TrialDurations)+posttrialtime+timeBinSize/2;


            %% Loading data from kilosort/phy easily
            myKsDir = fullfile(KilosortDir,MiceOpt{midx});
            subksdirs = dir(fullfile(myKsDir,thisdate,'**','Probe*')); %This changed because now I suddenly had 2 probes per recording
            if length(subksdirs)<1
                clear subksdirs
                subksdirs.folder = myKsDir; %Should be a struct array
                subksdirs.name = 'Probe0';
            end
            ProbeOpt = (unique({subksdirs(:).name}));
            Alignflag = 0; %only needs to be done once
            for probeid = 1:length(ProbeOpt)
                %Saving directory
                thisprobe = ProbeOpt{probeid}
                myKsDir = fullfile(KilosortDir,MiceOpt{midx},thisdate,thisprobe);

                % Check for multiple subfolders?
                subsesopt = dir(fullfile(myKsDir,'**','channel_positions.npy'));
                subsesopt = unique({subsesopt(:).folder});
                if isempty(subsesopt)
                    disp(['No data found in ' myKsDir ', continue...'])
                    continue
                end
                %Saving directory
                tmpfile = dir(fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,'MPEPData.mat'));
                if ~Redo && ~isempty(tmpfile) && datetime(tmpfile.date)>= FromDate
                    disp([MiceOpt{midx} ' ' thisdate   ' ' thisses ' already processed on ' tmpfile.date '... skipping'])
                    if CopyToXFileSpecific
                        if ~isdir(fullfile('\\znas.cortexlab.net\Lab\Share\Enny\MPEPData_Preprocessed',xfile,MiceOpt{midx},thisdate,thisses,thisprobe))
                            mkdir(fullfile('\\znas.cortexlab.net\Lab\Share\Enny\MPEPData_Preprocessed',xfile,MiceOpt{midx},thisdate,thisses,thisprobe))
                        end
                        copyfile(fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,'MPEPData.mat'),fullfile('\\znas.cortexlab.net\Lab\Share\Enny\MPEPData_Preprocessed',xfile,MiceOpt{midx},thisdate,thisses,thisprobe,'MPEPData.mat'))
                    end
                    continue
                end

                if ~Redo && exist(fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,'MPEPData.mat'))
                    if ~RedoAfterClustering || exist(fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,thisprobe,'CuratedResults.mat'))
                        continue
                    elseif RedoAfterClustering
                        myKsDir = fullfile(KilosortDir,MiceOpt{midx},thisdate);
                        myClusFile = dir(fullfile(myKsDir,'cluster_info.tsv'));
                        if isempty(myClusFile)
                            disp('This data is not yet curated with phy!!')
                            continue
                        end
                        NewHistologyNeeded = 1; %Automatically to 1 after RedoAfterClustering
                    end
                end

                %% Align to Trial Onset times
                if ~Alignflag
                    [starttrialidx,endtrialidx,trialid] = FindTrialOnAndOffSetsPhotoDiode(Actualtime,Timeline(:,ismember(AllInputs,'photoDiode')),TrialDurations);

                    TrialDurations(isnan(starttrialidx))=[];
                    endtrialidx(isnan(starttrialidx))=[];
                    starttrialidx(isnan(starttrialidx))=[];
                    if isempty(starttrialidx)
                        continue
                    end
                    if any((Actualtime(endtrialidx)-Actualtime(starttrialidx)) - TrialDurations'>1/tmSR*10)
                        disp(['flips slightly drifting... Average of ' num2str(nanmean((Actualtime(endtrialidx)-Actualtime(starttrialidx)) - TrialDurations')) 'sec'])

                    end
                    if ntrials ~= length(starttrialidx)
                        warning('Can''t find enough trials')
                        continue
                    end
                    Alignflag = 1;
                end
                %Saving directory
                if ~isfolder(fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe))
                    mkdir(fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe))
                end

                % Remove current processed data
                delete(fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,'*'))

                %                 %% Computing some useful details about spikes/neurons (like depths)
                myLFDir = fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},thisdate,'ephys');
                lfpD = dir(fullfile(myLFDir,'*','*','*.ap.*bin')); % ap file from spikeGLX specifically
                if isempty(lfpD)
                    disp('No data found.. skip')
                    continue
                end
                if length(lfpD)~=length(subksdirs)
                    if ~isempty(dir(fullfile(myLFDir,'*VR*','*','*.ap.*bin')))
                        lfpD = dir(fullfile(myLFDir,'*VR*','*','*.ap.*bin')); % ap file from spikeGLX specifically
                    elseif length(lfpD)~=length(subksdirs)
                        disp('Should be a different amount of probes?')
                    end
                end
                lfpD = lfpD(probeid);

                % Get information from meta file
                [Imecmeta] = ReadMeta2(lfpD.folder,'ap');
                lfpFs = str2num(Imecmeta.imSampRate);
                nChansInFile = strsplit(Imecmeta.acqApLfSy,',');  % neuropixels phase3a, from spikeGLX
                nChansInFile = str2num(nChansInFile{1})+1; %add one for sync


                %% Get cluster information
                PipelineParams.thisdate = thisdate;
                PipelineParams.SaveDir = fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe);
                try
                    [clusinfo, sp, Params]  = LoadPreparedClusInfo(subsesopt,PipelineParams);
                catch ME
                    disp(ME)
                    PipelineParams = ExtractKilosortData(subsesopt,PipelineParams);
                    [clusinfo, sp, Params]  = LoadPreparedClusInfo(subsesopt,PipelineParams);
                end
                % This extracts the parameters within clusinfo and sp
                % struct for further analysis
                ExtractFields({sp,clusinfo})

                %% load synchronization data
                SyncKSDataToTimeline


                % In this case, take only relevant recording sesion
                Good_IDx = find(Good_ID & ismember(RecSesID',recordingsessionidx));
                nclus = length(Good_IDx);
                if nclus < 2
                    disp('Less than 2 good units, skip...')
                    continue
                end
                if isfield(Params,'UnitMatch') && Params.UnitMatch == 1
                    if ~isfield(clusinfo,'UniqueID')
                        disp('No UniqueID... skip...')
                        continue
                    end
                    cluster_idUsed = UniqueID;
                    clu_sp = UniqClu;

                    if all(isnan(UniqueID(Good_IDx)))
                        disp('No UniqueID... skip')
                        cluster_idUsed = cluster_id;
                        clu_sp = clu;
                    elseif any(isnan(UniqueID(Good_IDx)))
                        keyboard
                    end
                else
                    cluster_idUsed = cluster_id;
                    clu_sp = clu;
                end

                %if necessary
                if syncchanmissing % There were some sessions where the clock was outputted from imec and this signal was also written on flipper. Try to extract that, combined with neuronal data
                    try
                        syncchanmissingTrySyncAgain
                    catch ME
                        disp(ME)
                        abortthissession = 1;
                    end

                end
                if abortthissession
                    continue
                end
                if ~syncchanmissing && ~any(~isnan(spikeTimesCorrected))
                    warning('No Spikes in this session... continue')
                    continue
                end

                %% Get Histology output
                %                 if strcmp(ProbeType{midx},'2_4S')
                %                     thisdate = []; % There's no data for the chronic mice in front of histology.
                %                 end
                clear Depth2Area
                GetHistologyOutput      %I created an extra function to have one line of code in the different scripts
                if ~histoflag
                    disp([MiceOpt{midx} thisdate thisses thisprobe 'No histology data...'])

                end
                thisdate = Dates4Mouse{didx}; % Reassign thisdate

                %% Spikes per trial
                SpikeIDx = arrayfun(@(X) find(spikeTimesCorrected>=Actualtime(starttrialidx(X))-pretrialtime&spikeTimesCorrected<=Actualtime(endtrialidx(X))+posttrialtime),1:ntrials,'UniformOutput',0);
                SpikeTrialID = nan(1,length(spikeTimesCorrected));
                SpikeTrialTime = nan(1,length(spikeTimesCorrected));
                for trid = 1:ntrials
                    SpikeTrialID(SpikeIDx{trid})=trid;
                    SpikeTrialTime(SpikeIDx{trid}) = spikeTimesCorrected(SpikeIDx{trid})-Actualtime(starttrialidx(trid));
                end
                %remove NaNs
                spikeDepths(isnan(SpikeTrialID))=[];
                SpikeTrialTime(isnan(SpikeTrialID))=[];
                clu_sp(isnan(SpikeTrialID)) = [];
                spikeTimesCorrected(isnan(SpikeTrialID))=[];
                RecSes(isnan(SpikeTrialID))=[];
                spikeShank(isnan(SpikeTrialID))=[];
                SpikeTrialID(isnan(SpikeTrialID))=[];

                %% Spike rate / histogram
                SpikeRatePerTP = arrayfun(@(Y) arrayfun(@(X) histcounts(SpikeTrialTime(SpikeTrialID == X & clu_sp'== cluster_idUsed(Y) & spikeShank' == clusinfo.Shank(Y) & RecSes' == clusinfo.RecSesID(Y)),...
                    timeEdges),1:ntrials,'UniformOutput',0),Good_IDx,'UniformOutput',0);
                SpikeRatePerTP = cat(1,SpikeRatePerTP{:});
                SpikeRatePerTP = cat(1,SpikeRatePerTP{:});
                SpikeRatePerTP = reshape(SpikeRatePerTP,length(Good_IDx),ntrials,[]); %Reshape to nclus, ntrials, ntp
                SpikeRatePerTP=SpikeRatePerTP./timeBinSize; %in spikes/sec
                SpikeRatePerTP = permute(SpikeRatePerTP,[2,3,1]); % Convert to trial,time,unit
                SpikeRatePerTP(isnan(SpikeRatePerTP))=0;
                %                 SpikeRatePerTP(repmat(sum(SpikeRatePerTP==0,3)==nclus,[1,1,nclus]))=nan;%
                %                 Fill with nans instead of 0 when longer %Unsure why I did
                %                 this?
                %%
                if ~any(SpikeRatePerTP(:)>0)
                    disp(['No spikes for ' MiceOpt{midx} ' ' thisdate ' ' thisses ', skip..'])
                    continue
                end


                %% Initialize save data
                MpepInfo = table;
                if size(cluster_idUsed,1)==1
                    cluster_idUsed=cluster_idUsed';
                end
                if size(cluster_id,1)==1
                    cluster_id=cluster_id';
                end
                if size(UniqueID,1)==1
                    UniqueID=UniqueID';
                end
                MpepInfo.ClusIDUsed = cluster_idUsed(Good_IDx);
                MpepInfo.ClusID = cluster_id(Good_IDx);
                MpepInfo.UniqueID = UniqueID(Good_IDx);

                MpepInfo.depth = depth(Good_IDx);
                MpepInfo.RecSes = clusinfo.RecSesID(Good_IDx);
                MpepInfo.Shank = clusinfo.Shank(Good_IDx);
                MpepInfo.SpikesPerSec = permute(SpikeRatePerTP,[3,1,2]);

                [sorteddepth,sortidx] = sort(depth(Good_IDx));

                AllCondNames={};
                for condid=1:ncond
                    conditionvals = arrayfun(@(X) [' ' ParDefs{X} '=' num2str(Pars(X,condid))],parrel,'UniformOutput',0);
                    conditionvals = {[conditionvals{:}]};
                    AllCondNames = {AllCondNames{:} conditionvals{1}};
                end

                TrialInfo = table;
                TrialInfo.CondNames = arrayfun(@(X) AllCondNames{X},condindx,'UniformOutput',0)';
                TrialInfo.CondIndx = condindx';
                if histoflag
                    areacol = clusinfo.Color(Good_IDx,:);
                    areatmp = clusinfo.Area(Good_IDx);
                    areatmp = strrep(areatmp,'/','');
                end
                %% Plot per condition across variables
                %                 figure('name','Average PSTH')
                [condsorted,sortid] = sort(condindx);
                %                 imagesc(newtimevec,condsorted,nanmean(SpikeRatePerTP(sortid,:,:),3))
                %                 colormap hot
                %                 ylabel('Condition')
                %                 xlabel('Time (s)')
                %                 saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,['AveragePSTH.fig']))
                %                 saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,['AveragePSTH.bmp']))
                %
                figure('name','Average PSTH','units','normalized','outerposition',[0 0 1 1])
                for parid = 1:length(parrel)
                    [CondOpt,idx,cdx] = unique(Pars(parrel(parid),:));
                    TrialsPerCond = nan(ceil(size(SpikeRatePerTP,1)./length(CondOpt)),length(CondOpt));
                    % Find trial index per condition
                    for cid = 1:length(CondOpt)
                        Cond2Take = find(Pars(parrel(parid),:)==CondOpt(cid)); %These are the conditions to take
                        Cond2Take(ismember(Cond2Take,ControlCond))=[]; %Remove control condition
                        TrialsPerCond(1:length(find(ismember(condindx,Cond2Take))),cid) = find(ismember(condindx,Cond2Take)); %These are the trials to take
                    end
                    TrialsPerCond(sum(isnan(TrialsPerCond),2)==size(TrialsPerCond,2),:)=[];
                    TrialsPerCond(TrialsPerCond==0)=nan;
                    subplot(1,length(parrel),parid)
                    imagesc(newtimevec,[],nanmean(SpikeRatePerTP(TrialsPerCond(~isnan(TrialsPerCond)),:,:),3))
                    set(gca,'YTick',[size(TrialsPerCond,1)/2:size(TrialsPerCond,1):length(TrialsPerCond(:))-size(TrialsPerCond,1)/2],...
                        'YTickLabel',CondOpt)
                    colormap hot
                    ylabel('Condition')
                    xlabel('Time (s)')
                    title([ParFullnames{parrel(parid)}])


                end
                saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,['AveragePSTH.fig']))
                saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,['AveragePSTH.bmp']))

                for parid = 1:length(parrel)
                    [CondOpt,idx,cdx] = unique(Pars(parrel(parid),:));
                    TrialsPerCond = nan(ceil(size(SpikeRatePerTP,1)./length(CondOpt)),length(CondOpt));
                    AllCondNames = {};

                    % Find trial index per condition
                    for cid = 1:length(CondOpt)
                        Cond2Take = find(Pars(parrel(parid),:)==CondOpt(cid)); %These are the conditions to take
                        Cond2Take(ismember(Cond2Take,ControlCond))=[]; %Remove control condition
                        TrialsPerCond(1:length(find(ismember(condindx,Cond2Take))),cid) = find(ismember(condindx,Cond2Take)); %These are the trials to take
                        conditionvals = [ParDefs{parrel(parid)} '=' num2str(CondOpt(cid))];
                        AllCondNames = {AllCondNames{:} conditionvals};
                    end
                    TrialsPerCond(sum(isnan(TrialsPerCond),2)==size(TrialsPerCond,2),:)=[];
                    TrialsPerCond(TrialsPerCond==0)=nan;

                    eval(['TrialInfo.' ParDefs{parrel(parid)} '= Pars(parrel(parid),condindx)'';'])
                    cols = jet(length(CondOpt));

                    figure('name',['PSTH per condition ' ParFullnames{parrel(parid)}],'units','normalized','outerposition',[0 0 1 1])
                    for cid = 1:length(CondOpt)
                        subplot(ceil(sqrt(length(CondOpt))),round(sqrt(length(CondOpt))),cid)
                        tmp = squeeze(nanmean(SpikeRatePerTP(TrialsPerCond(~isnan(TrialsPerCond(:,cid)),cid),:,:),1));
                        h=shadedErrorBar(newtimevec,nanmean(tmp,2),nanstd(tmp,[],2)./sqrt(nclus-1),'lineProps',{'color',cols(cid,:),'LineWidth',1.5});
                        ylabel('Spikes/sec')
                        xlabel('Time (s)')
                        title(AllCondNames{cid})
                    end
                    saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,[ParDefs{parrel(parid)} ' AveragePSTHPerCondition.fig']))
                    saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,[ParDefs{parrel(parid)} ' AveragePSTHPerCondition.bmp']))
                end
                %% TF analysis
                try
                    fb=cwtfilterbank('SamplingFrequency',1/timeBinSize,'SignalLength',length(newtimevec),'FrequencyLimits',freqlims,'VoicesPerOctave',nv);
                    figure
                    freqz(fb)
                    psi = wavelets(fb);
                    F = centerFrequencies(fb);
                    % Spikes
                    tmpspks = nanmean(nanmean(SpikeRatePerTP,1),3);
                    tmpspks(isnan(tmpspks))=0;
                    [csf,f] = cwt(double(tmpspks),'FilterBank',fb);

                    % Spikes
                    tmpspks = nanmean(nanmean(SpikeRatePerTP,1),3);
                    tmpspks(isnan(tmpspks))=0;
                    [csf,f] = cwt(double(tmpspks),'FilterBank',fb);
                    % amplitude?
                    Ampl= abs(csf);

                    %Time frequency for different conditions?
                    figure('name','Average amplitude across trials','units','normalized','outerposition',[0 0 1 1]);
                    subplot(2,1,1)
                    h=imagesc(newtimevec,[],nanmean(Ampl,3));
                    xlabel('Time (s)')
                    ylabel('Frequency Indx')
                    title('Amplitude')
                    colormap hot
                    set(gca,'ydir','normal')

                    makepretty

                    subplot(2,1,2)
                    avgAmpl = nanmean(Ampl,2);
                    plot(f,avgAmpl);
                    xlabel('Frequency (Hz)')
                    ylabel('Max Amplitude')
                    makepretty
                catch ME
                    disp(ME)
                end

                %%
                try
                    for parid = 1:length(parrel)
                        [CondOpt,idx,cdx] = unique(Pars(parrel(parid),:));
                        TrialsPerCond = nan(ceil(size(SpikeRatePerTP,1)./length(CondOpt)),length(CondOpt));
                        AllCondNames = {};

                        % Find trial index per condition
                        for cid = 1:length(CondOpt)
                            Cond2Take = find(Pars(parrel(parid),:)==CondOpt(cid)); %These are the conditions to take
                            Cond2Take(ismember(Cond2Take,ControlCond))=[]; %Remove control condition
                            TrialsPerCond(1:length(find(ismember(condindx,Cond2Take))),cid) = find(ismember(condindx,Cond2Take)); %These are the trials to take
                            conditionvals = [ParDefs{parrel(parid)} '=' num2str(CondOpt(cid))];
                            AllCondNames = {AllCondNames{:} conditionvals};
                        end
                        TrialsPerCond(sum(isnan(TrialsPerCond),2)==size(TrialsPerCond,2),:)=[];
                        TrialsPerCond(TrialsPerCond==0)=nan;

                        figure('name',['TF per condition ' ParDefs{parrel(parid)}],'units','normalized','outerposition',[0 0 1 1])
                        for cid = 1:length(CondOpt)
                            subplot(ceil(sqrt(length(CondOpt))),round(sqrt(length(CondOpt))),cid)
                            % Spikes
                            tmpspks = nanmean(nanmean(SpikeRatePerTP(TrialsPerCond(~isnan(TrialsPerCond(:,cid)),cid),:,:),1),3);
                            tmpspks(isnan(tmpspks))=0;
                            [csf,f] = cwt(double(tmpspks),'FilterBank',fb);
                            Ampl= abs(csf);


                            h=imagesc(newtimevec,[],Ampl);
                            xlabel('Time (s)')
                            colormap hot
                            title(AllCondNames{cid})
                            set(gca,'YTick','')
                            hold on

                            yyaxis right
                            avgAmpl = nanmean(Ampl,2);
                            plot(avgAmpl,f,'b-');
                            ylabel('Frequency (Hz)')
                            set(gca,'YScale','log')
                            makepretty
                        end
                        saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,[ParDefs{parrel(parid)} ' TFperCondition.fig']))
                        saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,[ParDefs{parrel(parid)} ' TFperCondition.bmp']))
                    end


                    %% Analyze without plotting for all units
                    if any(ismember(RelName,{'tf','ModFreq'})) % Frequency modulation!
                        parvec = find(ismember(ParDefs,{'tf','ModFreq'}));

                        figure('name',['Freq modulation across depth - Rayleigh'],'units','normalized','outerposition',[0 0 1 1]);

                        for parid = 1:length(parvec)
                            freqs = unique(Pars(parvec(parid),:))/10;
                            if length(freqs)==1
                                continue
                            end
                            Rvec = nan(nclus,length(freqs));
                            Thetavec = nan(nclus,length(freqs));
                            RayleighP = nan(nclus,length(freqs));
                            for freqid=1:length(freqs)

                                trialidx = find(ismember(condindx,find(Pars(parvec(parid),:)==freqs(freqid)*10)));
                                circledur = 1./freqs(freqid);
                                angles = arrayfun(@(Y) arrayfun(@(X) (newtimevec(SpikeRatePerTP(X,:,Y)>0)*2*pi)./circledur,trialidx,'UniformOutput',0),1:nclus,'UniformOutput',0);
                                angles = cellfun(@(Y) [Y{:}],angles,'UniformOutput',0); %Concatenate across trials
                                clusid = ~cell2mat(cellfun(@isempty,angles,'UniformOutput',0));
                                RayleighP(clusid,freqid)=cell2mat(cellfun(@(Y) circ_rtest(Y'),angles(clusid),'UniformOutput',0)); % Rayleightest
                                Thetavec(clusid,freqid) = cell2mat(cellfun(@(Y) circ_mean(Y'),angles(clusid),'UniformOutput',0)); % circular mean
                                Rvec(clusid,freqid) = cell2mat(cellfun(@(Y) circ_r(Y'),angles(clusid),'UniformOutput',0)); % strength

                            end
                            %Save
                            eval(['MpepInfo.' ParDefs{parvec(parid)} '_Rvec=Rvec;']) % ModFreq should have 2 points per variable to fit in large mpep table
                            eval(['MpepInfo.' ParDefs{parvec(parid)} '_CircMean=Thetavec;'])
                            eval(['MpepInfo.' ParDefs{parvec(parid)} '_RayleighP=RayleighP;'])


                            subplot(1,5,[(parid-1)*2+1 parid*2])
                            % Plot this
                            Pvals = abs(log10(RayleighP));
                            imagesc(freqs,sorteddepth,Pvals(sortidx,:),[0 10])
                            xlabel('Frequency')
                            if parid==1
                                ylabel('Depth (micron from tip)')
                            else
                                set(gca,'YTickLabel',[])
                            end
                            set(gca,'YDir','normal','XTick',freqs)
                            title(ParDefs{parvec(parid)})
                            colormap hot
                            freezeColors
                        end
                        if histoflag
                            [areaopt,idx,udx] = unique(areatmp(sortidx),'stable');
                            areacolsort = areacol(sortidx);
                            subplot(1,5,5)
                            h=imagesc(udx);

                            colormap(cat(1,areacolsort{idx}))
                            set(gca,'YDir','normal')
                            set(gca,'YTick',idx,'YTickLabel',areaopt,'XTickLabel',[])

                        end

                        saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,['Freq modulation across depth - Rayleigh.fig']))
                        saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,['Freq modulation across depth - Rayleigh.bmp']))

                    end
                catch ME
                    disp(ME)
                end


                %% Plot per unit
                nexample = 5

                TotalSpikeRate = squeeze(nansum(nansum(SpikeRatePerTP,1),2));
                if nexample>sum(TotalSpikeRate>5000)
                    nexample=sum(TotalSpikeRate>5000);
                end
                exampleunits = datasample(find(TotalSpikeRate>5000),nexample,'replace',false);
                for uid = 1:nexample
                    unitid = cluster_idUsed(Good_IDx(exampleunits(uid)));
                    depthhere = depth(Good_IDx(exampleunits(uid)));
                    if histoflag
                        unitname = ['Unit ' num2str(unitid) ', depth=' num2str(depthhere) ', ' areatmp{exampleunits(uid)}];
                    else
                        unitname = ['Unit ' num2str(unitid) ', depth=' num2str(depthhere)];
                    end

                    for parid = 1:length(parrel)
                        [CondOpt,idx,cdx] = unique(Pars(parrel(parid),:));
                        TrialsPerCond = nan(ceil(size(SpikeRatePerTP,1)./length(CondOpt)),length(CondOpt));
                        AllCondNames = {};
                        cols = jet(length(CondOpt));
                        AllCols = nan(3,ceil(size(SpikeRatePerTP,1)./length(CondOpt)),length(CondOpt));
                        % Find trial index per condition
                        for cid = 1:length(CondOpt)
                            Cond2Take = find(Pars(parrel(parid),:)==CondOpt(cid)); %These are the conditions to take
                            Cond2Take(ismember(Cond2Take,ControlCond))=[]; %Remove control condition
                            TrialsPerCond(1:length(find(ismember(condindx,Cond2Take))),cid) = find(ismember(condindx,Cond2Take)); %These are the trials to take
                            conditionvals = [ParDefs{parrel(parid)} '=' num2str(CondOpt(cid))];
                            AllCondNames = {AllCondNames{:} conditionvals};
                            AllCols(:,1:length(find(ismember(condindx,Cond2Take))),cid) = repmat(cols(cid,:)',[1,length(find(ismember(condindx,Cond2Take)))]);
                        end
                        AllCols(:,sum(isnan(TrialsPerCond),2)==size(TrialsPerCond,2),:)=[];
                        AllCols = reshape(AllCols,3,[]);
                        TrialsPerCond(sum(isnan(TrialsPerCond),2)==size(TrialsPerCond,2),:)=[];
                        TrialsPerCond(TrialsPerCond==0)=nan;

                        AllCols = AllCols(:,~isnan(TrialsPerCond(:)));
                        TrialsPerCond = TrialsPerCond(~isnan(TrialsPerCond));
                        try
                            figure('name',[unitname ParDefs{parrel(parid)}],'units','normalized','outerposition',[0 0 1 1])
                            subplot(2,1,1)
                            hold on

                            arrayfun(@(X) scatter(newtimevec(SpikeRatePerTP(TrialsPerCond(X),:,exampleunits(uid))>0),repmat(X,[1,sum(SpikeRatePerTP(TrialsPerCond(X),:,exampleunits(uid))>0)]),8,AllCols(:,X)','filled'),1:length(TrialsPerCond),'UniformOutput',0)
                            xlabel('Time (s)')
                            ylabel('Trial (sorted by condition)')
                            line([0 0],get(gca,'ylim'),'color',[0 0 0],'LineStyle','--')

                            TrialsPerCond = reshape(TrialsPerCond,[],length(CondOpt));
                            subplot(2,1,2)
                            hold on
                            clear h
                            for cid=1:length(CondOpt)
                                tmp = SpikeRatePerTP(TrialsPerCond(:,cid),:,exampleunits(uid));
                                h(cid)=shadedErrorBar(newtimevec,smooth(nanmean(tmp,1)+max(get(gca,'ylim'))*0.9),smooth(nanstd(tmp,[],1)./sqrt(size(tmp,1)-1)),'lineProps',{'color',cols(cid,:),'LineWidth',1.5});
                            end
                            xlabel('time (s)')
                            ylabel('Spks/Sec')
                            set(gca,'YTick','')
                            line([0 0],get(gca,'ylim'),'color',[0 0 0],'LineStyle','--')

                            legend([h.mainLine],AllCondNames)
                            saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,[ParDefs{parrel(parid)} unitname '_PSTH.fig']))
                            saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,[ParDefs{parrel(parid)} unitname '_PSTH.bmp']))
                        catch ME
                            disp(ME)
                        end
                    end
                    if any(ismember(RelName,{'tf','ModFreq'})) % Frequency modulation!
                        parvec = find(ismember(ParDefs,{'tf','ModFreq'}));
                        for parid = 1:length(parvec)
                            freqs = unique(Pars(parvec(parid),:))/10;
                            if freqs==0
                                continue
                            end

                            figure('name',[unitname ParDefs{parvec(parid)} ' Frequency Modulation'],'units','normalized','outerposition',[0 0 1 1])
                            for freqid=1:length(freqs)

                                trialidx = find(ismember(condindx,find(Pars(parvec(parid),:)==freqs(freqid)*10)));
                                circledur = 1./freqs(freqid);
                                angles = arrayfun(@(X) (newtimevec(SpikeRatePerTP(X,:,exampleunits(uid))>0)*2*pi)./circledur,trialidx,'UniformOutput',0);
                                if isempty([angles{:}])
                                    continue
                                end

                                subplot(length(freqs),2,(2*freqid)-1)
                                try
                                    histogram([angles{:}]',[0:0.2*pi:max([angles{:}])])
                                    ylims = get(gca,'ylim');
                                    arrayfun(@(X) patch([X X+pi X+pi X],[min(ylims) min(ylims) max(ylims) max(ylims)],[0 1 0],'FaceAlpha',0.2,'EdgeColor','none'),0:2*pi:max([angles{:}]),'UniformOutput',0)
                                    set(gca,'XTickLabel',cellfun(@(X) num2str(round(str2num(X)./(2*pi).*circledur*100)./100),(get(gca,'XTickLabel')),'UniformOutput',0))
                                    xlabel('Time (s)')
                                    hold on
                                    ylabel(['nr spikes at phase ' ParDefs{parvec(parid)}])
                                    title([num2str(freqs(freqid)) 'Hz'])
                                    makepretty

                                    subplot(length(freqs),2,2*freqid)
                                    circ_plot([angles{:}]','hist',[],20,true,true,'linewidth',2,'color','r')
                                    p=circ_rtest([angles{:}]');
                                    title(['p=' num2str(p)])
                                    makepretty
                                catch ME
                                    disp(ME)
                                end
                            end

                            saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,[unitname ParDefs{parvec(parid)} ' Frequency Modulation.fig']))
                            saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,[unitname ParDefs{parvec(parid)} ' Frequency Modulation.bmp']))

                        end
                    end
                end

                %% Analyze  for different depths

                if 0
                    dptstp = 500; %in um
                    depthvec = min(spikeDepths):dptstp:max(spikeDepths);
                    CondOpt=unique(condindx);


                    SpikeRatePerDepth = arrayfun(@(Y) arrayfun(@(X) histcounts(SpikeTrialTime((SpikeTrialID == X)' & spikeDepths>=depthvec(Y)-dptstp/2&spikeDepths<=depthvec(Y)+dptstp/2),...
                        timeEdges),1:ntrials,'UniformOutput',0),1:length(depthvec),'UniformOutput',0);
                    SpikeRatePerDepth = cat(1,SpikeRatePerDepth{:});
                    SpikeRatePerDepth = cat(1,SpikeRatePerDepth{:});
                    SpikeRatePerDepth = reshape(SpikeRatePerDepth,length(newtimevec),length(depthvec),[]);
                    storedepth = [];
                    cols = summer(length(CondOpt));
                    for parid = 1:length(parrel)
                        [CondOpt,idx,cdx] = unique(Pars(parrel(parid),:));
                        TrialsPerCond = nan(ceil(size(SpikeRatePerTP,1)./length(CondOpt)),length(CondOpt));
                        AllCondNames = {};
                        cols = jet(length(CondOpt));
                        AllCols = nan(3,ceil(size(SpikeRatePerTP,1)./length(CondOpt)),length(CondOpt));
                        % Find trial index per condition
                        for cid = 1:length(CondOpt)
                            Cond2Take = find(Pars(parrel(parid),:)==CondOpt(cid)); %These are the conditions to take
                            Cond2Take(ismember(Cond2Take,ControlCond))=[]; %Remove control condition
                            TrialsPerCond(1:length(find(ismember(condindx,Cond2Take))),cid) = find(ismember(condindx,Cond2Take)); %These are the trials to take
                            conditionvals = [ParDefs{parrel(parid)} '=' num2str(CondOpt(cid))];
                            AllCondNames = {AllCondNames{:} conditionvals};
                            AllCols(:,1:length(find(ismember(condindx,Cond2Take))),cid) = repmat(cols(cid,:)',[1,length(find(ismember(condindx,Cond2Take)))]);
                        end
                        AllCols(:,sum(isnan(TrialsPerCond),2)==size(TrialsPerCond,2),:)=[];
                        AllCols = reshape(AllCols,3,[]);
                        TrialsPerCond(sum(isnan(TrialsPerCond),2)==size(TrialsPerCond,2),:)=[];
                        TrialsPerCond = reshape(TrialsPerCond,[],length(CondOpt));



                        for dptid=1:length(depthvec)
                            largefig = figure('name',['Spikes depth ' num2str(depthvec(dptid)) 'um' ParDefs{parrel(parid)}],'units','normalized','outerposition',[0 0 1 1]);

                            % Extract spike indices for this depth and
                            % condition
                            idx = arrayfun(@(X) ismember(SpikeTrialID,X)' & spikeDepths>=depthvec(dptid)-dptstp/2&spikeDepths<=depthvec(dptid)+dptstp/2,TrialsPerCond(:),'UniformOutput',0);
                            % Extract spike times
                            spktms = (cellfun(@(X) SpikeTrialTime(X),idx,'UniformOutput',0));
                            subplot(2,1,1)
                            hold on
                            arrayfun(@(X) scatter(spktms{X},depthvec(dptid)+X-1,10,AllCols(:,X)','filled'),1:length(spktms))
                            xlabel('Time (s)')
                            ylabel('Depth and trials (sorted by condition)')
                            line([0 0],get(gca,'ylim'),'color',[0 0 0],'LineStyle','--')

                            subplot(2,1,2)
                            hold on
                            clear h
                            for cid=1:length(CondOpt)
                                tmp = squeeze(SpikeRatePerDepth(:,dptid,TrialsPerCond(:,cid)));
                                h(cid)=shadedErrorBar(newtimevec,smooth(nanmean(tmp,2)+depthvec(dptid)./25+cid),smooth(nanstd(tmp,[],2)./sqrt(size(tmp,2)-1)),'lineProps',{'color',cols(cid,:),'LineWidth',1.5});
                            end
                            drawnow
                            xlabel('time (s)')
                            ylabel('Spks/Sec')
                            set(gca,'YTick','')
                            line([0 0],get(gca,'ylim'),'color',[0 0 0],'LineStyle','--')

                            legend([h.mainLine],AllCondNames)

                            saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,[ParDefs{parrel(parid)} num2str(depthvec(dptid)) 'um' '_PSTH.fig']))
                            saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,[ParDefs{parrel(parid)} num2str(depthvec(dptid)) 'um' '_PSTH.bmp']))
                            close(largefig)
                        end

                    end
                end
                clear SpikeRatePerDepth

                %% Orientation tuning?
                if any(ismember(RelName,{'ori'}))
                    parvec = find(ismember(ParDefs,{'ori'}));

                    figure('name',['Orientation tuning across depth'],'units','normalized','outerposition',[0 0 1 1]);

                    for parid = 1:length(parvec)
                        oris = unique(Pars(parvec(parid),:));

                        SpikesPerOri = nan(nclus,ceil(size(SpikeRatePerTP,1)./length(oris)),length(oris));
                        for Oriid=1:length(oris)
                            %orientation, contrast not 0!
                            trialidx = find(ismember(condindx,find(Pars(parvec(parid),:)==oris(Oriid) & nanmean(Pars(ismember(ParDefs,{'cr','cg','cb'}),:),1)~=0)));
                            SpikesPerOri(:,1:length(trialidx),Oriid) = squeeze(nanmean(SpikeRatePerTP(trialidx,newtimevec>0&newtimevec<=nanmin(TrialDurations),:),2))'; %Average over stimulus presentation
                        end

                        % OSI
                        %Prefered orientation - odd trials
                        [~,PrefOriIdx] = (arrayfun(@(X) nanmax(nanmean(SpikesPerOri(X,1:2:end,:),2),[],3),1:nclus,'UniformOutput',0));
                        PrefOriIdx = cell2mat(PrefOriIdx);
                        PrefOri = oris(PrefOriIdx);
                        Plus90Ori = mod(PrefOri+90,360);
                        Plus90OriIdx = cell2mat(arrayfun(@(X) find(oris==X),Plus90Ori,'UniformOutput',0));

                        OSI = cell2mat(arrayfun(@(X) (nanmean(SpikesPerOri(X,2:2:end,PrefOriIdx(X)),2)-nanmean(SpikesPerOri(X,2:2:end,Plus90OriIdx(X)),2))./...
                            (nanmean(SpikesPerOri(X,2:2:end,PrefOriIdx(X)),2)+nanmean(SpikesPerOri(X,2:2:end,Plus90OriIdx(X)),2)),1:nclus,'UniformOutput',0));
                        % DSI
                        Plus180Ori = mod(PrefOri+180,360);
                        Plus180OriIdx = cell2mat(arrayfun(@(X) find(oris==X),Plus180Ori,'UniformOutput',0));
                        DSI = cell2mat(arrayfun(@(X) (nanmean(SpikesPerOri(X,2:2:end,PrefOriIdx(X)),2)-nanmean(SpikesPerOri(X,2:2:end,Plus180OriIdx(X)),2))./...
                            (nanmean(SpikesPerOri(X,2:2:end,PrefOriIdx(X)),2)+nanmean(SpikesPerOri(X,2:2:end,Plus180OriIdx(X)),2)),1:nclus,'UniformOutput',0));

                        figure('units','normalized','outerposition',[0 0 1 1],'units','normalized','outerposition',[0 0 1 1]); subplot(2,2,2)
                        histogram(OSI)
                        makepretty
                        title('OSI')
                        subplot(2,2,4)
                        histogram(DSI)
                        makepretty
                        title('DSI')
                        OSI = OSI';
                        DSI = DSI';

                        % Significant Rayleigh?

                        spks = squeeze(round(nansum(SpikesPerOri(:,:,:),2)));
                        angles = arrayfun(@(Y) arrayfun(@(X) repmat(oris(X),[1,spks(Y,X)]),1:length(oris),'UniformOutput',0),...
                            1:nclus,'UniformOutput',0);
                        angles = cat(1,angles{:});
                        unitidx = find(sum(cell2mat(cellfun(@isempty,angles,'UniformOutput',0)),2)<length(oris));
                        pRayleigh = nan(1,nclus);
                        pRayleigh(unitidx)=cell2mat(arrayfun(@(Y) circ_rtest([angles{Y,:}]'./360*2*pi),unitidx,'UniformOutput',0));

                        subplot(2,2,[1,3])
                        cols = repmat([0 0 0],nclus,1);
                        cols(pRayleigh<0.05,:)=repmat([1 0 0],sum(pRayleigh<0.05),1);
                        scatter(abs(log(pRayleigh)),depth(Good_IDx),10,cols,'filled')
                        hold on
                        pRayleigh = pRayleigh';
                        PrefOri = PrefOri';
                        %Save
                        eval(['MpepInfo.' ParDefs{parvec(parid)} '_SpikesPerOri=SpikesPerOri;'])
                        eval(['MpepInfo.' ParDefs{parvec(parid)} '_OSI=OSI;'])
                        eval(['MpepInfo.' ParDefs{parvec(parid)} '_DSI=DSI;'])
                        eval(['MpepInfo.' ParDefs{parvec(parid)} '_pRayleigh=pRayleigh;'])
                        eval(['MpepInfo.' ParDefs{parvec(parid)} '_PrefOri=PrefOri;'])


                    end


                    saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,['OrientationModulation.fig']))
                    saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,['OrientationModulation.bmp']))

                    % Few example cells
                    for uid = 1:nexample
                        unitid = cluster_idUsed(Good_IDx(exampleunits(uid)));
                        depthhere = depth(Good_IDx(exampleunits(uid)));
                        if histoflag
                            unitname = ['Unit ' num2str(unitid) ', depth=' num2str(depthhere) ', ' areatmp{exampleunits(uid)}];
                        else
                            unitname = ['Unit ' num2str(unitid) ', depth=' num2str(depthhere)];
                        end

                        figure('name',unitname,'units','normalized','outerposition',[0 0 1 1])

                        subplot(2,1,1)
                        shadedErrorBar(oris,squeeze(nanmean(SpikesPerOri(exampleunits(uid),:,:),2)),squeeze(nanstd(SpikesPerOri(exampleunits(uid),:,:),[],2)))
                        xlabel('Orientation')
                        ylabel('Spks/Sec')
                        title(['OSI=' num2str(OSI(exampleunits(uid))) ', DSI=' num2str(DSI(exampleunits(uid)))])
                        makepretty

                        % Polar Plot
                        subplot(2,1,2)
                        spks = round(nansum(squeeze(SpikesPerOri(exampleunits(uid),:,:)),1));
                        angles = arrayfun(@(X) repmat(oris(X),[1,spks(X)]),1:length(oris),'UniformOutput',0);
                        angles = [angles{:}];
                        if ~isempty(angles)
                            circ_plot(angles'./360*2*pi,'hist',[],30,true,true,'linewidth',2,'color','r')
                            p=circ_rtest(angles'./360*2*pi);
                            title(['p=' num2str(p)])
                        end
                        saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,[unitname '_PSTH.fig']))
                        saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,[unitname '_PSTH.bmp']))



                    end

                end

                %% Frequency tuning?
                if any(ismember(RelName,{'SoundFreq'}))
                    parvec = find(ismember(ParDefs,{'SoundFreq'}));

                    figure('name',['Frequency tuning across depth'],'units','normalized','outerposition',[0 0 1 1]);

                    for parid = 1:length(parvec)
                        oris = unique(Pars(parvec(parid),:));

                        SpikesPerOri = nan(nclus,ceil(size(SpikeRatePerTP,1)./length(oris)),length(oris));
                        for Oriid=1:length(oris)
                            %Frequency, Amplitude not 0!
                            trialidx = find(ismember(condindx,find(Pars(parvec(parid),:)==oris(Oriid) & nanmean(Pars(ismember(ParDefs,{'ampl'}),:),1)~=0)));
                            SpikesPerOri(:,1:length(trialidx),Oriid) = squeeze(nanmean(SpikeRatePerTP(trialidx,newtimevec>0&newtimevec<=nanmin(TrialDurations),:),2))'; %Average over stimulus presentation
                        end

                        % FMI
                        %Prefered frequency - odd trials
                        [~,PrefFreqIdx] = (arrayfun(@(X) nanmax(nanmean(SpikesPerOri(X,1:2:end,:),2),[],3),1:nclus,'UniformOutput',0));
                        PrefFreqIdx = cell2mat(PrefFreqIdx);
                        PrefFreq = oris(PrefFreqIdx);

                        [~,NPrefFreqIdx] = (arrayfun(@(X) nanmin(nanmean(SpikesPerOri(X,1:2:end,:),2),[],3),1:nclus,'UniformOutput',0));
                        NPrefFreqIdx = cell2mat(NPrefFreqIdx);
                        NPrefFreq = oris(NPrefFreqIdx);

                        FMI = cell2mat(arrayfun(@(X) (nanmean(SpikesPerOri(X,2:2:end,PrefFreqIdx(X)),2)-nanmean(SpikesPerOri(X,2:2:end,NPrefFreqIdx(X)),2))./...
                            (nanmean(SpikesPerOri(X,2:2:end,PrefFreqIdx(X)),2)+nanmean(SpikesPerOri(X,2:2:end,NPrefFreqIdx(X)),2)),1:nclus,'UniformOutput',0));

                        histogram(FMI)
                        makepretty
                        title('Frequency Modulation Index')

                        FMI = FMI';
                        PrefFreq = PrefFreq';
                        NPrefFreq = NPrefFreq';

                        %Save
                        eval(['MpepInfo.' ParDefs{parvec(parid)} '_SpikesPerFreq=SpikesPerOri;'])
                        eval(['MpepInfo.' ParDefs{parvec(parid)} '_FMI=FMI;'])
                        eval(['MpepInfo.' ParDefs{parvec(parid)} '_PrefFreq=PrefFreq;'])
                        eval(['MpepInfo.' ParDefs{parvec(parid)} '_NPrefFreq=NPrefFreq;'])


                    end


                    saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,['FrequencyModulation.fig']))
                    saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,['FrequencyModulation.bmp']))

                    % Few example cells
                    figure('name',unitname,'units','normalized','outerposition',[0 0 1 1])
                    for uid = 1:nexample
                        unitid = cluster_idUsed(Good_IDx(exampleunits(uid)));
                        depthhere = depth(Good_IDx(exampleunits(uid)));
                        if histoflag
                            unitname = ['Unit ' num2str(unitid) ', depth=' num2str(depthhere) ', ' areatmp{exampleunits(uid)}];
                        else
                            unitname = ['Unit ' num2str(unitid) ', depth=' num2str(depthhere)];
                        end

                        subplot(nexample,1,uid)
                        shadedErrorBar(oris,squeeze(nanmean(SpikesPerOri(exampleunits(uid),:,:),2)),squeeze(nanstd(SpikesPerOri(exampleunits(uid),:,:),[],2)))
                        xlabel('frequency')
                        ylabel('Spks/Sec')
                        title(['FMI=' num2str(FMI(exampleunits(uid)))])
                        makepretty


                    end

                end

                %% Optogenetic effect?
                if any(ismember(RelName,{'amp'})) && nclus>10


                    % why are there some trials with more activity?
                    OptoStart = unique(Pars(strcmp(ParDefs,'tstart'),:)./1000);
                    OptoStop = OptoStart+unique(Pars(strcmp(ParDefs,'tend'),:)./1000);
                    % Unit/sorting management
                    inclUnits = find(sum(squeeze(nanmean(SpikeRatePerTP,2))==0,1)<0.5*ntrials); %Exclude units without much activity
                    [depthsorted,sortidx] = sort(depth(Good_IDx(inclUnits)));
                    % 
                    % % Z-score
                    % basez = squeeze(nanmean(SpikeRatePerTP(:,:,inclUnits),1));
                    % stdz = squeeze(nanstd(SpikeRatePerTP(:,:,inclUnits),[],1));
                    % zSc = (permute(SpikeRatePerTP(:,:,inclUnits),[2,3,1])-repmat(basez,[1,1,ntrials]))./repmat(stdz,[1,1,ntrials]); % Z-score per trial and cluster over time
                    % zSc = permute(zSc,[3,1,2]);

                    basez = squeeze(nanmean(reshape(SpikeRatePerTP(:,newtimevec<OptoStart|newtimevec>OptoStop+0.5,inclUnits),[],length(inclUnits)),1));
                    stdz = squeeze(nanstd(reshape(SpikeRatePerTP(:,newtimevec<OptoStart|newtimevec>OptoStop+0.5,inclUnits),[],length(inclUnits)),[],1));
                    zSc = (permute(SpikeRatePerTP(:,:,inclUnits),[3,1,2])-basez')./stdz'; % Z-score per trial and cluster over time
                    zSc = permute(zSc,[2,3,1]);

                    
                    figure('name','General spiking','units','normalized','outerposition',[0 0 1 1])
                    subplot(2,2,1)
                    imagesc(squeeze(nanmean(zSc(:,:,sortidx),2))',[-.5 .5])
                    colormap redblue
                    xlabel('Trial')
                    ylabel('Unit (sorted by depth)')
                    title('Z-scored activit')
                    makepretty

                    subplot(2,2,2)
                    imagesc([],newtimevec,squeeze(nanmean(zSc(:,:,sortidx),3))',[-.5 .5])
                    colormap redblue
                    xlabel('Trial')
                    ylabel('Time (s)')
                    title('Z-scored activit')
                    makepretty

                    subplot(2,2,3)
                    imagesc(newtimevec,[],squeeze(nanmean(zSc(:,:,sortidx),1)),[-0.5 0.5])
                    hold on
                    line([OptoStart OptoStart],get(gca,'ylim'),'color',[0 0 0],'LineStyle','--')
                    line([OptoStop OptoStop],get(gca,'ylim'),'color',[0 0 0],'LineStyle','--')
                    colormap redblue
                    ylabel('Unit (sorted by depth)')
                    xlabel('Time (s)')
                    title('Z-scored activit')
                    makepretty

                    subplot(2,2,4)
                    imagesc([],[],(squeeze(nanmean(SpikeRatePerTP(:,newtimevec>OptoStart&newtimevec<OptoStop,inclUnits),2))'-squeeze(nanmean(nanmean(SpikeRatePerTP(sortid,newtimevec<OptoStart,inclUnits),2),1))),[-40 40])
                    colormap redblue
                    title('Spks/s')
                    xlabel('Trial')
                    ylabel('Time (s)')
                    makepretty

                    % Any units with modulation?
                    %                     tmp = reshape(SpikeRatePerTP,[],nclus)';
                    %                     zSc = reshape(zSc',size(SpikeRatePerTP));


                    amps = unique(Pars(strcmp(ParDefs,'amp'),:));
                    shapevec = Pars(strcmp(ParDefs,'shape'),:);
                    shapes = unique(shapevec);
                    ShapeOpt = {'Ramp','Step','Half-sine','Tapered'};
                    ShapeOpt = ShapeOpt(unique(shapes));
                    if numel(shapes)==1 % Switched at some point to using same shape, different pattern of activity
                        shapevec = Pars(strcmp(ParDefs,'freq'),:);
                        shapes = unique(shapevec);
                        ShapeOpt = arrayfun(@(X) ['freq=' num2str(X/10)],shapes,'Uni',0);
                    end            
                    cols = cat(3,jet(length(amps)),jet(length(amps)).*0.8);

                    SpikesPerAmp = nan(length(inclUnits),ceil(size(SpikeRatePerTP,1)./length(amps)),length(amps),numel(shapes));
                    AmplitudeVec = nan(1,ntrials);
                    ShapeVec = nan(1,ntrials);
                    AllCols = nan(3,ntrials);
                    for shapeid = 1:numel(shapes)
                        for Oriid=1:length(amps)
                            %orientation, contrast not 0!
                            trialidx = find(ismember(condindx,find(Pars(strcmp(ParDefs,'amp'),:)==amps(Oriid) & shapevec==shapes(shapeid))));
                            AmplitudeVec(trialidx) = amps(Oriid);
                            ShapeVec(trialidx) = shapes(shapeid);
                            AllCols(:,trialidx) = repmat(squeeze(cols(Oriid,:,shapeid))',1,length(trialidx));
                            SpikesPerAmp(:,1:length(trialidx),Oriid,shapeid) = squeeze(nanmean(zSc(trialidx,newtimevec>OptoStart&newtimevec<=OptoStop,:),2))'; %Average over stimulus presentation
                        end
                    end
                    [Ampsorted,sortid] = sort(AmplitudeVec,'ascend');

                    %                     BaseMod = squeeze(nanmean(zSc(:,newtimevec<OptoStart,:),2));

                    if histoflag
                        colperunit = clusinfo.Color(Good_IDx(inclUnits),:);
                    else
                        colperunit = copper(length(inclUnits));
                    end
                    colperunit = colperunit(sortidx,:);

                    % Correlation of OptoMod & AmplitudeVec
                    % Any units with opto modulation??
                    OptoMod = squeeze(nanmean(zSc(:,newtimevec>OptoStart&newtimevec<OptoStop,:),2));
                    [tmpcor, tmppval] = corr(OptoMod,AmplitudeVec');

                    % Average OptoMod per amplitude
                    OptoModMean = arrayfun(@(X) squeeze(nanmean(nanmean(zSc(AmplitudeVec==X,newtimevec>OptoStart&newtimevec<OptoStop,:),2),1)),amps,'Uni',0);
                    OptoModMean = cat(2,OptoModMean{:});

                    OptoModStd = arrayfun(@(X) squeeze(nanstd(nanmean(zSc(AmplitudeVec==X,newtimevec>OptoStart&newtimevec<OptoStop,:),2),[],1))./sqrt(sum(AmplitudeVec==X)-1),amps,'Uni',0);
                    OptoModStd = cat(2,OptoModStd{:});


                    figure('name','ResponseDuringOpto','units','normalized','outerposition',[0 0 1 1])
                    % Neurons modulating up
                    subplot(1,4,1)
                    inclUnits = find(tmppval<0.05 & tmpcor>0);
                    hold on
                    for id = 1:length(inclUnits)
                        errorbar(amps,OptoModMean(inclUnits(id),:)',OptoModStd(inclUnits(id),:)','color',colperunit(inclUnits(id),:))
                    end
                    ylim([-0.4 0.4])
                    xlabel('Light Intensity (a.u.)')
                    ylabel('z-scored activity')
                    title('Positively correlated')
                    makepretty

                    subplot(1,4,2)
                    inclUnits = find(tmppval<0.05 & tmpcor<0);
                    hold on
                    for id = 1:length(inclUnits)
                        errorbar(amps,OptoModMean(inclUnits(id),:)',OptoModStd(inclUnits(id),:)','color',colperunit(inclUnits(id),:))
                    end
                    ylim([-0.4 0.4])
                    xlabel('Light Intensity (a.u.)')
                    ylabel('z-scored response')
                    title('Negatively correlated')
                    makepretty

                    subplot(1,4,3)
                    inclUnits = find(tmppval>=0.05);
                    hold on
                    for id = 1:length(inclUnits)
                        errorbar(amps,OptoModMean(inclUnits(id),:)',OptoModStd(inclUnits(id),:)','color',colperunit(inclUnits(id),:))
                    end
                    ylim([-0.4 0.4])
                    xlabel('Light Intensity (a.u.)')
                    ylabel('z-scored response')
                    title('Not correlated')
                    makepretty

                    % refind all units 2 iclude
                    inclUnits = find(sum(squeeze(nanmean(SpikeRatePerTP,2))==0,1)<0.5*ntrials); %Exclude units without much activity

                    subplot(1,4,4)
                    h=imagesc(1,depthsorted,[1:length(inclUnits)]')
                    colormap(copper(length(inclUnits)))
                    ylabel('Depth')
                    set(gca,'ydir','normal','XTick',[])

                    saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,['AmplitudeEffectErrorBar.fig']))
                    saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,['AmplitudeEffectErrorBar.bmp']))

                    figure('units','normalized','outerposition',[0 0 1 1]);
                    h1 = subplot(1,2,1);
                    imagesc(amps,depthsorted,OptoModMean(sortidx,:),[-0.3 0.3])
                    xlabel('Amplitude (a.u.)')
                    ylabel('Depth (micron)')
                    colormap redblue
                    makepretty
                    colorbar
                    title('Z-scored activity during opto ON')
                    set(gca,'ydir','normal')

                    h2 = subplot(1,2,2);
                    p = polyfit(depth(Good_IDx(inclUnits)),tmpcor,1);
                    y = polyval(p,depth(Good_IDx(inclUnits)));

                    scatter(tmpcor,depth(Good_IDx(inclUnits)),14,[0 0 0],'filled')
                    hold on
                    plot(y,depth(Good_IDx(inclUnits)),'r--')
                    ylabel('Depth (um)')
                    xlabel('correlation with Amplitude')

                    linkaxes([h1 h2],'y')


                    %                     % Example Units
                    %                     try
                    %                         exampleunits = randsample(find(abs(tmpcor)>0.5),1);
                    %                     catch
                    %                         try
                    %                         exampleunits = randsample(find(~isnan(tmpcor)),1);
                    %                         catch
                    %                             exampleunits = randsample(inclUnits,1);
                    %                         end
                    %                     end
                    %                     exampleunits = inclUnits(exampleunits); % IN bigger dataset
                    %
                    %                     for uid = 1:length(exampleunits)
                    %                         if histoflag
                    %                             unitname = ['Unit ' num2str(exampleunits(uid)) ', depth=' num2str(depth(Good_IDx(exampleunits(uid)))) ', ' areatmp{exampleunits(uid)}];
                    %                         else
                    %                             unitname = ['Unit ' num2str(exampleunits(uid)) ', depth=' num2str(depth(Good_IDx(exampleunits(uid)))) ];
                    %                         end
                    %
                    %                         subplot(2,2,4)
                    %                         hold on
                    %                         arrayfun(@(X) scatter(newtimevec(SpikeRatePerTP(sortid(X),:,exampleunits(uid))>0),repmat(X,[1,sum(SpikeRatePerTP(sortid(X),:,exampleunits(uid))>0)]),4,AllCols(:,sortid(X))','filled'),1:length(sortid),'UniformOutput',0)
                    %                         line([OptoStart OptoStart],get(gca,'ylim'),'color',[1 0 0],'LineStyle','--')
                    %                         line([OptoStop OptoStop],get(gca,'ylim'),'color',[1 0 0],'LineStyle','--')
                    %                         xlabel('Time (s)')
                    %                         ylabel('Trial (sorted by condition)')
                    % %
                    % %                         for Oriid=1:length(amps)
                    % %                             trialidx = find(ismember(condindx,find(Pars(strcmp(ParDefs,'amp'),:)==amps(Oriid))));
                    % %                             h=plot(newtimevec,nanmean(SpikeRatePerTP(trialidx,:,exampleUnits(uid)),1)+2.5*(Oriid-1));
                    % %                             hold on
                    % %                         end
                    % %                         xlabel('Time (s)')
                    % %                         ylabel('Spikes/sec')
                    %                         title(unitname)
                    %                         makepretty
                    %
                    %                     end
                    saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,['AmplitudeEffect.fig']))
                    saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,['AmplitudeEffect.bmp']))




                    %
                    try
                        exampleunits = randsample(find(abs(tmpcor)>0.3),5);
                    catch
                        try
                            exampleunits = randsample(find(~isnan(tmpcor)),5);
                        catch
                            exampleunits = randsample(inclUnits,6);
                        end
                    end
                    exampleunits = inclUnits(exampleunits); % IN bigger dataset

                    figure('name','Examples','units','normalized','outerposition',[0 0 1 1])
                    for uid = 1:length(exampleunits)
                        if histoflag
                            unitname = ['Unit ' num2str(clusinfo.cluster_id(Good_IDx(exampleunits(uid)))) ', depth=' num2str(depth(Good_IDx(exampleunits(uid)))) ', ' areatmp{exampleunits(uid)}];
                        else
                            unitname = ['Unit ' num2str(clusinfo.cluster_id(Good_IDx(exampleunits(uid)))) ', depth=' num2str(depth(Good_IDx(exampleunits(uid)))) ];
                        end

                        subplot(2,5,uid)
                        hold on
                        arrayfun(@(X) scatter(newtimevec(SpikeRatePerTP(sortid(X),:,exampleunits(uid))>0),repmat(X,[1,sum(SpikeRatePerTP(sortid(X),:,exampleunits(uid))>0)]),4,AllCols(:,sortid(X))','filled'),1:length(sortid),'UniformOutput',0)
                        xlabel('Time (s)')
                        ylabel('Trial (sorted by condition)')
                        line([OptoStart OptoStart],get(gca,'ylim'),'color',[1 0 0],'LineStyle','--')
                        line([OptoStop OptoStop],get(gca,'ylim'),'color',[1 0 0],'LineStyle','--')
                        title(unitname)
                        makepretty

                        subplot(2,5,uid+5)
                        hold on
                        clear h
                        for cid=1:length(amps)
                            tmp = squeeze(SpikeRatePerTP(AmplitudeVec == amps(cid),:,exampleunits(uid)));
                            h(cid)=shadedErrorBar(newtimevec,smooth(squeeze(nanmean(tmp,1)),2)+cid*50,smooth(squeeze(nanstd(tmp,[],1))./sqrt(size(tmp,2)-1)),'lineProps',{'color',cols(cid,:,1),'LineWidth',1.5});
                        end
                        drawnow
                        xlabel('time (s)')
                        ylabel('Spks/Sec')
                        set(gca,'YTick','')
                        line([OptoStart OptoStart],get(gca,'ylim'),'color',[1 0 0],'LineStyle','--')
                        line([OptoStop OptoStop],get(gca,'ylim'),'color',[1 0 0],'LineStyle','--')

                    end
                    legend([h.mainLine],arrayfun(@(X) ['Amp = ' num2str(X./1000)],amps,'Uni',0))

                    saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,['Amplitude_PSTH.fig']))
                    saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,['Amplitude_PSTH.bmp']))

                    %% Fit a proper effect?
                    % Average OptoMod per amplitude
                    OptoModMean = arrayfun(@(X) squeeze(nanmean(nanmean(zSc(AmplitudeVec==X,newtimevec>OptoStart&newtimevec<OptoStop,:),2),1)),amps,'Uni',0);
                    OptoModMean = cat(2,OptoModMean{:});

                    OptoModStd = arrayfun(@(X) squeeze(nanstd(nanmean(zSc(AmplitudeVec==X,newtimevec>OptoStart&newtimevec<OptoStop,:),2),[],1))./sqrt(sum(AmplitudeVec==X)-1),amps,'Uni',0);
                    OptoModStd = cat(2,OptoModStd{:});

                    expFun = @(p,a) p(1)+p(2).*(a./(p(3)+a));%+p(3); % For spatial decay
                    opts = optimset('Display','off');
                    GoodnessofFit = nan(1,numel(inclUnits));
                    Params = nan(3,numel(inclUnits));


                    for uid = 1:numel(inclUnits)
                        % tmpfig = figure;

                        try
                            p = lsqcurvefit(expFun,[nanmin(OptoModMean(uid,:)) 0 median(amps)],double(amps),OptoModMean(uid,:),[-inf -10 0],[inf 10 max(amps)],opts);
                            GoodnessofFit(uid) = nansum((OptoModMean(uid,:)-expFun(p,double(amps)))./max(abs(OptoModMean(uid,:))).^2);
                            Params(:,uid) = p;

                            % 
                            % shadedErrorBar(amps,OptoModMean(uid,:),OptoModStd(uid,:))
                            % hold on
                            % plot(double(amps),expFun(p,double(amps)),'r-')
                            % drawnow
                            % pause(0.2)

                        catch ME
                            disp(ME)
                        end
                        % close(tmpfig)

                    end
                    figure('name','Parameter fits')
                    subplot(2,2,1)
                    histogram(Params(1,:))
                    title('Baseline response')

                    subplot(2,2,2)
                    histogram(Params(2,:))
                    title('K')


                    subplot(2,2,3)
                    histogram(Params(3,:))
                    title('A50')

                    subplot(2,2,4)
                    histogram(GoodnessofFit)
                    title('Goodness of fit')
                    saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,['ParameterFits.fig']))
                    saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,['ParameterFits.bmp']))


                    figure('name','Scatter vs Kappa')
                    subplot(2,2,1)
                    scatter(Params(2,:),tmpcor)
                    xlabel('Kappa')
                    ylabel('Correlation')
                    subplot(2,2,2)
                    scatter(Params(2,:),Params(3,:))
                    xlabel('Kappa')
                    ylabel('A50')
                    subplot(2,2,3)
                    scatter(Params(1,:),Params(2,:))
                    xlabel('R0')
                    ylabel('Kappa')
                    subplot(2,2,4)
                    scatter(Params(1,:),Params(3,:))
                    ylabel('A50')
                    xlabel('R0')

                    saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,['ScattervsKappa.fig']))
                    saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,['ScattervsKappa.bmp']))

                    %
                    figure('name','Effect on activity')
                    for shid = 1:numel(shapes)
                        subplot(2,2,(shid-1)*2+1)
                        tmp = reshape(SpikesPerAmp(tmppval<0.05&tmpcor>0,:,:,shid),sum(tmppval<0.05&tmpcor>0),[]);
                        InclHere = sum(~isnan(tmp),1)~=0;
                        ampshere = repmat(amps,size(SpikesPerAmp,2),1);
                        ampshere = ampshere(:);
                        imagesc(ampshere(InclHere),[],tmp(:,InclHere),[-1 1]);
                        colormap redblue
                        xlabel('Amplitude')
                        title([ShapeOpt{shid} ' activated by light'])'
                        colorbar
                        makepretty
                        offsetAxes

                        subplot(2,2,(shid-1)*2+2)
                        tmp = reshape(SpikesPerAmp(tmppval<0.05&tmpcor<0,:,:,shid),sum(tmppval<0.05&tmpcor<0),[]);
                        InclHere = sum(~isnan(tmp),1)~=0;
                        ampshere = repmat(amps,size(SpikesPerAmp,2),1);
                        ampshere = ampshere(:);
                        imagesc(ampshere(InclHere),[],tmp(:,InclHere),[-1 1]);
                        colormap redblue
                        xlabel('Amplitude')
                        title([ShapeOpt{shid} ' suppressed by light'])'
                        colorbar
                        makepretty
                        offsetAxes
                    end
                    saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,['AffectedPopulation.fig']))
                    saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,['AffectedPopulation.bmp']))

                    %% RSP neurons specifically?
                    if histoflag
                        RSPIdx = find(contains(areatmp(inclUnits),'rsp'));
                        notRSPIdx = find(~contains(areatmp(inclUnits),'rsp'));
                        clusIds = clusinfo.cluster_id(Good_IDx(inclUnits(RSPIdx)));

                        wdur = tdfread(fullfile(myKsDir,'cluster_waveform_dur.tsv'));
                        wdur = wdur.waveform_dur(ismember(wdur.cluster_id,clusIds));

                        FSNeuron = wdur<=400; % Putative fast spiking neuron (GABA)

                        figure('name','RSP opto versus neuron type')

                        subplot(1,2,1)
                        scatter(wdur,Params(2,RSPIdx),25,[1 0 0].*FSNeuron,'filled')
                        hold on
                        line(get(gca,'xlim'),[0 0],'color',[0.2 0.2 0.2])
                        ylabel('Kappa')
                        xlabel('waveform duration')
                        makepretty
                        offsetAxes

                        subplot(1,2,2)
                        scatter(wdur,tmpcor(RSPIdx),25,[1 0 0].*FSNeuron,'filled')
                        hold on
                        line(get(gca,'xlim'),[0 0],'color',[0.2 0.2 0.2])

                        xlabel('waveform duration')
                        ylabel('Correlation with amplitude')
                        makepretty
                        offsetAxes
                        saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,['WaveformVsKapaRSPOnly.fig']))
                        saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,['WaveformVsKapaRSPOnly.bmp']))




                        %
                        figure('name','Effect on RSP and no RSP activity')
                        for shid = 1:2
                            subplot(2,2,(shid-1)*2+1)
                            tmp = reshape(SpikesPerAmp(RSPIdx,:,:,shid),numel(RSPIdx),[]);
                            InclHere = sum(~isnan(tmp),1)~=0;
                            ampshere = repmat(amps,size(SpikesPerAmp,2),1);
                            ampshere = ampshere(:);
                            imagesc(ampshere(InclHere),[],tmp(:,InclHere),[-1 1]);
                            colormap redblue
                            xlabel('Amplitude')
                            title([ShapeOpt{shid} ' RSP neurons'])'
                            colorbar
                            makepretty
                            offsetAxes

                            subplot(2,2,(shid-1)*2+2)
                            tmp = reshape(SpikesPerAmp(notRSPIdx,:,:,shid),numel(notRSPIdx),[]);
                            InclHere = sum(~isnan(tmp),1)~=0;
                            ampshere = repmat(amps,size(SpikesPerAmp,2),1);
                            ampshere = ampshere(:);
                            imagesc(ampshere(InclHere),[],tmp(:,InclHere),[-1 1]);
                            colormap redblue
                            xlabel('Amplitude')
                            title([ShapeOpt{shid} ' Not RSP neurons'])'
                            colorbar
                            makepretty
                            offsetAxes
                        end
                        saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,['RSPvsNotRSP.fig']))
                        saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,['RSPvsNotRSP.bmp']))






                    end


                    tmpcortable = nan(nclus,1);
                    tmpcortable(inclUnits) = tmpcor;
                    MpepInfo.OptoCorrelation = tmpcortable;
                    tmpcortable = nan(nclus,1);
                    tmpcortable(inclUnits) = tmppval;
                    MpepInfo.OptoCorrelationPval = tmpcortable;
                    tmpcortable = nan(nclus,1);
                    tmpcortable(inclUnits) = Params(1,:);
                    MpepInfo.R0 = tmpcortable;
                    tmpcortable = nan(nclus,1);
                    tmpcortable(inclUnits) = Params(2,:);
                    MpepInfo.Kappa = tmpcortable;
                    tmpcortable = nan(nclus,1);
                    tmpcortable(inclUnits) = Params(3,:);
                    MpepInfo.A50 = tmpcortable;


                end



                %% Save data
                if exist('Depth2AreaPerUnit')
                    save(fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,'MPEPData.mat'),'MpepInfo','TrialInfo','Protocol','Depth2AreaPerUnit','clusinfo','-v7.3')
                else
                    save(fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,'MPEPData.mat'),'MpepInfo','TrialInfo','Protocol','clusinfo','-v7.3')
                end

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
                clear clu_sp
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
                clear UniqueID

                close all
            end
        end
    end
end
return

%% table across mice
PipelineParams=OriSetting;
if RedoTable || ~exist(fullfile(SaveDir,'AllMiceMPEPData.mat'))
    createtable = 1;
    sesstypes={};
    AllTrialInfo = table;
    for midx = 1:length(MiceOpt)
        Dates4Mouse = DateOpt{midx};
        if PipelineParams.UnitMatch %% USE UNIQUE_ID
            tmp = load(fullfile(SaveDir,MiceOpt{midx},'UnitMatch','UnitMatch.mat'));
            MatchTable = tmp.MatchTable;
            UMParam = tmp.UMparam;
            AllKSDir = UMParam.KSDir;
            UniqueIDConversion = tmp.UniqueIDConversion;
            RecSesID = UniqueIDConversion.recsesAll;
            OriID =  UniqueIDConversion.OriginalClusID;
            UniqueID =  UniqueIDConversion.UniqueID;
        else
            keyboard
            adduid=0;
        end
        for didx = 1:length(Dates4Mouse)
            % Within folders, look for 'RF mapping sessions'
            thisdate = Dates4Mouse{didx};

            subsess = dir(fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},Dates4Mouse{didx}));
            subsess(1:2) = []; %remove '.' and '..'
            flag = 0;
            mpepsess = [];
            for sesidx=1:length(subsess)
                if strcmp(subsess(sesidx).name,'leftovers')
                    continue
                end
                listfiles = dir(fullfile(subsess(sesidx).folder,subsess(sesidx).name,'*.p')); %Find protocol file
                if any(~cell2mat(cellfun(@isempty,cellfun(@(X) strfind(X,'.p'),{listfiles(:).name},'UniformOutput',0),'UniformOutput',0)))
                    idx = find(~cell2mat(cellfun(@isempty,cellfun(@(X) strfind(X,'.p'),{listfiles(:).name},'UniformOutput',0),'UniformOutput',0)));
                    % read which x-file was used
                    if length(idx)>1
                        disp('Too many files?!')
                        keyboard
                    end
                    fileID = fopen(fullfile(listfiles(idx).folder,listfiles(idx).name));
                    A = fscanf(fileID,'%c');
                    fclose(fileID);

                    if ~any(strfind(A,'SparseNoise'))
                        SesName = strsplit(A,'.x');
                        SesName = SesName{1};
                        mpepsess = [mpepsess {subsess(sesidx).name}];
                        sesstypes = {sesstypes{:} SesName};
                        flag = 1;
                        continue
                    else
                        disp('RF session, use different script... continue')
                    end
                else
                    continue
                end
            end

            if ~flag
                continue
            end

            for sesidx = 1:length(mpepsess)
                thisses = mpepsess{sesidx};
                allfiles = dir(fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,'*','MPEPData.mat'));
                if isempty(allfiles)
                    continue
                end
                for pib=1:length(allfiles)
                    % Load VR Data
                    tmp = load(fullfile(allfiles(pib).folder,allfiles(pib).name));
                    ProbeId= strsplit(allfiles(pib).folder,'Probe');
                    TrialInfo = tmp.TrialInfo;
                    ntrials = length(TrialInfo.CondIndx);
                    TrialInfo.Mouse = repmat(MiceOpt{midx},ntrials,1);
                    TrialInfo.Date = repmat(thisdate,ntrials,1);
                    TrialInfo.Session = repmat(thisses,ntrials,1);
                    thisprobe = ['Probe' ProbeId{2}];
                    TrialInfo.Probe = repmat(['Probe' ProbeId{2}],ntrials,1);

                    % Add area to mpepinfo table
                    try
                        tmp.MpepInfo.Area = arrayfun(@(X) tmp.Depth2AreaPerUnit.Area{ismember(tmp.Depth2AreaPerUnit.cluster_idUsed,tmp.MpepInfo.ClusID(X)) & ismember(tmp.Depth2AreaPerUnit.Shank,tmp.MpepInfo.Shank(X))},1:length(tmp.MpepInfo.ClusID),'Uni',0)';
                        tmp.MpepInfo.Color = arrayfun(@(X) tmp.Depth2AreaPerUnit.Color{ismember(tmp.Depth2AreaPerUnit.cluster_idUsed,tmp.MpepInfo.ClusID(X)) & ismember(tmp.Depth2AreaPerUnit.Shank,tmp.MpepInfo.Shank(X))},1:length(tmp.MpepInfo.ClusID),'Uni',0)';

                    catch
                        tmp.MpepInfo.Area = arrayfun(@(X) tmp.Depth2AreaPerUnit.Area{ismember(tmp.Depth2AreaPerUnit.cluster_idUsed,tmp.MpepInfo.ClusID(X)) & ismember(tmp.Depth2AreaPerUnit.Depth,tmp.MpepInfo.depth(X))},1:length(tmp.MpepInfo.ClusID),'Uni',0)';
                        tmp.MpepInfo.Color = arrayfun(@(X) tmp.Depth2AreaPerUnit.Color{ismember(tmp.Depth2AreaPerUnit.cluster_idUsed,tmp.MpepInfo.ClusID(X)) & ismember(tmp.Depth2AreaPerUnit.Depth,tmp.MpepInfo.depth(X))},1:length(tmp.MpepInfo.ClusID),'Uni',0)';

                    end
                    tmp.MpepInfo.Area  = strrep(tmp.MpepInfo.Area,'/','');

                    if any(~cellfun(@length, tmp.MpepInfo.Color))
                        keyboard
                    end
                    % Add session&mouse&probe information
                    tmp.MpepInfo.Mouse = repmat(MiceOpt{midx},length(tmp.MpepInfo.ClusID),1);
                    tmp.MpepInfo.Date = repmat(thisdate,length(tmp.MpepInfo.ClusID),1);
                    tmp.MpepInfo.Session = repmat(thisses,length(tmp.MpepInfo.ClusID),1);
                    tmp.MpepInfo.Probe = repmat(thisprobe,length(tmp.MpepInfo.ClusID),1);

                    % Add UniqueID
                    if PipelineParams.UnitMatch %% USE UNIQUE_ID from UnitMatch
                        rechere = find(cellfun(@(X) contains(X,thisdate) & contains(X,thisprobe),AllKSDir));
                        if length(rechere)>1  % Select correct session
                            removeidx = zeros(1,length(rechere));
                            for id = 1:length(rechere)
                                tmpfile = matfile(fullfile(AllKSDir{rechere(id)},'PreparedData.mat'));
                                clusinfo = tmpfile.clusinfo;
                                if sum(clusinfo.Good_ID)~=length(tmp.RFInfo.ClusID)
                                    removeidx(id)=1;
                                end
                            end
                            rechere(removeidx)=[];
                        end
                        tmp.MpepInfo.UniqueID = arrayfun(@(X) UniqueID(RecSesID==rechere & OriID==X),tmp.MpepInfo.ClusID);
                    else
                        keyboard
                        tmp.MpepInfo.UniqueID = (1:length(tmp.MpepInfo.ClusID))+adduid;
                        if pib==1
                            adduid = adduid+length(tmp.MpepInfo.ClusID);
                        end
                    end
                    % Exclude spks per sec here
                    tmp.MpepInfo.SpikesPerSec = [];
                    tmpvariablenames = tmp.MpepInfo.Properties.VariableNames;

                    % Save in big struct
                    if createtable
                        createtable = 0;
                        AllMPEPDat = tmp.MpepInfo;
                        %                         AllMPEPDat = tall(AllMPEPDat); %Make tall to prevent memory issues
                        AllTrialInfo = TrialInfo;
                        AllVariables = AllMPEPDat.Properties.VariableNames;
                    else
                        includedvariables = ismember(tmpvariablenames,AllVariables);
                        if any(~includedvariables)
                            idxvec = find(~includedvariables);
                            clear tmptbl
                            tmptbl = table;
                            tmptbl.ClusID = gather(AllMPEPDat.ClusID);
                            for idx = 1:length(idxvec)
                                % new variables, add to table
                                eval(['tmpvar = tmp.MpepInfo.' tmpvariablenames{idxvec(idx)} ';'])
                                eval(['tmptbl.' tmpvariablenames{idxvec(idx)} ' = nan(gather(size(AllMPEPDat,1)),size(tmpvar,2));'])
                            end

                            % Add here
                            AllMPEPDat = innerjoin(AllMPEPDat,tall(tmptbl),'LeftKeys','ClusID','RightKeys','ClusID');
                            AllVariables = AllMPEPDat.Properties.VariableNames;
                        end

                        % create tmptable that matches the size of
                        % AllMPEPDat
                        nonincludedvariables = find(~ismember(AllVariables,tmpvariablenames));
                        tmptbl = tmp.MpepInfo;
                        for idx = 1:length(nonincludedvariables)
                            % new variables, add to table
                            sz = size(tmptbl);
                            eval(['tmpvar = AllMPEPDat.' AllVariables{nonincludedvariables(idx)} ';'])
                            eval(['tmptbl.' AllVariables{nonincludedvariables(idx)} ' = nan(size(tmptbl,1),gather(size(tmpvar,2)));'])
                        end
                        tmpvariablenames = tmptbl.Properties.VariableNames;


                        VariableIdx = cell2mat(cellfun(@(X) find(ismember(tmpvariablenames,X)),AllVariables,'UniformOutput',0));
                        tmptbl = tmptbl(:,VariableIdx);% In right order
                        % Concatenate
                        try
                            %                             GeTblClasses = varfun(@class,gather(AllMPEPDat),'OutputFormat','cell');
                            %                             scintTblClasses = varfun(@class,tmptbl,'OutputFormat','cell');
                            AllMPEPDat = cat(1,AllMPEPDat,tmptbl);

                            %                             AllMPEPDat = cat(1,AllMPEPDat,tall(tmptbl));
                            %                         AllTrialInfo = cat(1,AllTrialInfo,TrialInfo);
                        catch ME
                            disp(ME)
                            keyboard
                        end


                    end
                end
            end
        end
    end
    save(fullfile(SaveRFDir,'AllMiceMPEPData.mat'),'AllMPEPDat','-v7.3')
else
    load(fullfile(SaveRFDir,'AllMiceMPEPData.mat'))
end
%% Get information on areas in dataset
areatmp = gather(AllMPEPDat.Area);
[areaopt,aidx,audx] = unique(areatmp,'stable');
uniquearean = length(areaopt);
% Use Allen Brain Atlas to find all areas in the recordings
% Add units to area specific cell
% areatmp = areatmp(~cellfun(@isempty,areatmp));
atlastable = readtable(fullfile(AllenCCFPath,'structure_tree_safe_2017.csv'));
AutomaticAREASOfInterest = {};
AutomaticAREASOfInterestFullName = {};
automaticareasofinterestid = nan(1,uniquearean); %Index for which area
newareaname = cell(1,length(areatmp));
newareaabrev = cell(1,length(areatmp));
% Make names in table similar to those used here
atlastable.acronym = lower(atlastable.acronym);
atlastable.acronym= strrep(atlastable.acronym,'/','');
for areaid=1:length(areaopt)
    %Find structure_id_path of area
    structure_id_path = atlastable.structure_id_path{find(ismember(lower(atlastable.acronym),areaopt{areaid}))};
    % Cut off last part
    parts = strsplit(structure_id_path,'/');
    %     parts(cellfun(@isempty,parts))=[];
    if length(parts)<=4
        tmpareaname = areaopt{areaid}; %if we cannot go further up
    else
        newstructure_id_path = fullfile(parts{2:end-2});
        newstructure_id_path = ['/' strrep(newstructure_id_path,'\','/') '/'];
        % Find area with this new structure id
        tmpareaname = atlastable.acronym{find(ismember(atlastable.structure_id_path,newstructure_id_path))};
    end
    %Does this one already exist in Automatic AREAS Of Interest?
    if ~any(ismember(AutomaticAREASOfInterest,tmpareaname))
        %no, then create
        AutomaticAREASOfInterest = {AutomaticAREASOfInterest{:} tmpareaname};
        AutomaticAREASOfInterestFullName = {AutomaticAREASOfInterestFullName{:} atlastable.name{find(ismember(atlastable.structure_id_path,newstructure_id_path))}};
    end
    % Index correctly
    automaticareasofinterestid(areaid) = find(ismember(AutomaticAREASOfInterest,tmpareaname));
    newareaabrev(ismember(areatmp,areaopt{areaid})) = {tmpareaname};
    newareaname(ismember(areatmp,areaopt{areaid}))={atlastable.name{find(ismember(atlastable.structure_id_path,newstructure_id_path))}}; %Save out per unit
end
if any(ismember(AutomaticAREASOfInterest,'root'))
    AutomaticAREASOfInterestFullName{ismember(AutomaticAREASOfInterest,'root')}='root';
end


%% Summary
VisModPval = gather(AllMPEPDat.tf_RayleighP);
AudiModPval = gather(AllMPEPDat.ModFreq_RayleighP);

% Significance:
VisModPvalLog = abs(log10(VisModPval));
% VisModPvalLog(VisModPvalLog>=10)=10;
AudiModPvalLog = abs(log10(AudiModPval));
% AudiModPvalLog(AudiModPvalLog>=10)=10;
areanumber = length(AutomaticAREASOfInterestFullName);
color = cellfun(@(X) hex2rgb(X),gather(AllMPEPDat.Color),'UniformOutput',0);
color = cat(1,color{:});

ColPerArea = colorcube(areanumber+4);
ColPerArea = ColPerArea(2:end-3,:);

%%
figure('name','VisualVsAuditorymodulation','units','normalized','outerposition',[0 0 1 1])
clear h
include = true(1,areanumber);

pvalperunit = nan(areanumber,round(length(newareaname)*0.8),2); %areaXunitXvis/audit
for areaid=1:areanumber
    if strcmp(AutomaticAREASOfInterestFullName{areaid},'root')
        include(areaid)=0;
        continue
    end

    h(areaid) = scatter(nanmean(VisModPvalLog(ismember(newareaname,AutomaticAREASOfInterestFullName{areaid}),:),2),nanmean(AudiModPvalLog(ismember(newareaname,AutomaticAREASOfInterestFullName{areaid}),:),2),...
        15,ColPerArea(areaid,:),'filled');
    tmp = nanmean(VisModPvalLog(ismember(newareaname,AutomaticAREASOfInterestFullName{areaid}),:),2);
    pvalperunit(areaid,1:length(tmp),1)=tmp;
    tmp = nanmean(AudiModPvalLog(ismember(newareaname,AutomaticAREASOfInterestFullName{areaid}),:),2);
    pvalperunit(areaid,1:length(tmp),2)=tmp;

    hold on
end
lim = max([max(get(gca,'ylim')) max(get(gca,'xlim'))]);
xlim([0 lim])
ylim([0 lim])
line([0 lim],[0 lim],'color',[0 0 0],'LineWidth',1.5)
axis square
xlabel('Visual Modulation')
ylabel('Auditory Modulation')
title('log-p values')
legend([h(include)],AutomaticAREASOfInterestFullName(include))

%% violin
figure('name','Pvalue Distribution','units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1)
boxplot(pvalperunit(include,:,1)','plotstyle','compact','colors',ColPerArea(include,:))
makepretty
ylabel('Log Pvalue')
set(gca,'XTickLabel',[])
title('Visual')


subplot(2,1,2)
boxplot(pvalperunit(include,:,2)','plotstyle','compact','colors',ColPerArea(include,:),'labels',AutomaticAREASOfInterest(include));
makepretty
ylabel('Log Pvalue')
title('Auditory')

%% Number and percentage of units
thresh = abs(log10(0.01));
nsig = nan(areanumber,2,2,2); %areaXtotal/sigXvis/audit
include = true(1,areanumber);

for areaid=1:areanumber
    if strcmp(AutomaticAREASOfInterestFullName{areaid},'root')
        include(areaid)=0;
        continue
    end
    tmp1 = VisModPvalLog(ismember(newareaname,AutomaticAREASOfInterestFullName{areaid}),:);
    nsig(areaid,2,1,:) = size(tmp1,1);
    nsig(areaid,1,1,:) = sum(tmp1>thresh,1);
    tmp2 = AudiModPvalLog(ismember(newareaname,AutomaticAREASOfInterestFullName{areaid}),:);
    nsig(areaid,2,2,:) = size(tmp2,1);
    nsig(areaid,1,2,:) = sum(tmp2>thresh,1);
end
figure('name','Significant Modulation','units','normalized','outerposition',[0 0 1 1])
for id = 1:2
    ax1(id) = subplot(2,2,id);
    bar(nsig(include,:,1,id),'stacked')
    title(['Visual ' num2str(id) ' Modulation '])
    makepretty
    set(gca,'XTickLabel',[])
    ylabel('Nr Units')
    ax2(id) = subplot(2,2,id+2);
    bar(nsig(include,:,2,id),'stacked')
    title(['Auditory ' num2str(id) ' Modulation '])
    makepretty
    if id==2
        legend({'Significant','Total'})
    end
    set(gca,'XTick',1:sum(include),'XTickLabel',AutomaticAREASOfInterest(include),'XTickLabelRotation',90)
    ylabel('Nr Units')
    linkaxes([ax1 ax2])
end


%% Tuning consistency for UniqueID
UniqueID = gather(AllMPEPDat.UniqueID);
Date = gather(AllMPEPDat.Date);
Session = gather(AllMPEPDat.Session);
Mouse = gather(AllMPEPDat.Mouse);
%%
figure('name','FrequencyTuning consistency','units','normalized','outerposition',[0 0 1 1])

for midx = 1:length(MiceOpt)
    % All rows for this mouse
    rowidx = find(ismember(cellstr(Mouse),MiceOpt{midx}));
    %
    %         RayleighP(clusid,freqid) % Rayleightest
    %         Thetavec(clusid,freqid) % circular mean
    %         Rvec(clusid,freqid)  % strength

    UniqueUnitsOpt = unique(UniqueID(rowidx));
    UniqueDateOpt = unique(cellstr(Date(rowidx,:)));
    UniqueSesOpt = unique(cellstr(Session(rowidx,:)));
    subplot(ceil(sqrt(length(MiceOpt))),round(sqrt(length(MiceOpt))),midx)
    hold on
    cols = color(rowidx,:);
    %     tf_Rvec = gather(AllMPEPDat.tf_Rvec(rowidx,:));
    RayleighP = gather(AllMPEPDat.tf_RayleighP(rowidx,:));
    for uid=1:length(UniqueUnitsOpt)
        unitIdx=find(UniqueID(rowidx)==UniqueUnitsOpt(uid));
        tmpval = RayleighP(unitIdx,:);
        dateidx = cellfun(@(X) find(ismember(UniqueDateOpt,X)),cellstr(Date(rowidx(unitIdx),:)));
        sesidx = cellfun(@(X) find(ismember(UniqueSesOpt,X)),cellstr(Session(rowidx(unitIdx),:)));

        plot(tmpval(:,1),tmpval(:,2),'-','color',cols(uid,:))
    end
    xlabel('RayleighP_Stimulus1')
    ylabel('RayleighP_stimulus2')
end
%% Example scatter
if 0
    nexample = 5;
    areasofinterest = {'ca','rspagl'};
    audPThresh = abs(log10(0.5));
    visPThresh = abs(log10(0.5));
    idx = find(ismember(newareaabrev,areasofinterest)' & nanmean(VisModPvalLog,2)>visPThresh & nanmean(AudiModPvalLog,2)>audPThresh);
    exampleunits = datasample(idx,nexample,'replace',false);

    % Information from cluster table
    clusterid = gather(AllMPEPDat.ClusID);
    SpikesPerSec = gather(AllMPEPDat.SpikesPerSec);
    Mouse = gather(AllMPEPDat.Mouse);
    Date = gather(AllMPEPDat.Date);
    Session = gather(AllMPEPDat.Session);
    Probe = gather(AllMPEPDat.Probe);

    % Information from trial table
    cid = AllTrialInfo.CondIndx;
    CondNames = AllTrialInfo.CondNames;
    [condopt,idx,udx] = unique(cid);
    CondNames = CondNames(idx);
    cols = jet(length(CondNames));

    freqs = [3.3 6.1 5.0 7.1];
    for uid = 1:nexample
        % Grab trial information
        thismouse = Mouse(exampleunits(uid),:);
        thisdate = Date(exampleunits(uid),:);
        thisses = Session(exampleunits(uid),:);

        condindx = AllTrialInfo.CondIndx(ismember(cellstr(AllTrialInfo.Mouse),thismouse)& ismember(cellstr(AllTrialInfo.Date),thisdate) & ismember(cellstr(AllTrialInfo.Session),thisses) & ismember(cellstr(AllTrialInfo.Probe),'Probe0'));
        ntrials = length(condindx);
        [condsorted,sortid] = sort(condindx);

        figure('name',[thismouse ' ' thisdate ' ' thisses ': Unit ' num2str(clusterid(exampleunits(uid))) areatmp{exampleunits(uid)}],'units','normalized','outerposition',[0 0 1 1])

        tmpspks = squeeze(SpikesPerSec(exampleunits(uid),:,:));
        subplot(2,1,1)
        hold on
        arrayfun(@(X) scatter(newtimevec(tmpspks(sortid(X),:)>0),repmat(X,[1,sum(tmpspks(sortid(X),:)>0)]),8,cols(condsorted(X),:),'filled'),1:ntrials,'UniformOutput',0)
        %     arrayfun(@(X) scatter(newtimevec(tmpspks(sortid(X),:),repmat(X,[1,sum(tmpspks(sortid(X),:))]),8,cols(condsorted(X),:),'filled'),1:ntrials,'UniformOutput',0)
        xlabel('Time (s)')
        ylabel('Trials')
        line([0 0],get(gca,'ylim'),'color',[0 0 0],'LineStyle','--')

        subplot(2,1,2)
        hold on
        clear h
        for cid=1:ncond
            tmp = tmpspks(condindx==cid,:);
            h(cid)=shadedErrorBar(newtimevec,smooth(nanmean(tmp,1)+max(get(gca,'ylim'))*0.9),smooth(nanstd(tmp,[],1)./sqrt(size(tmp,1)-1)),'lineProps',{'color',cols(cid,:),'LineWidth',1.5});
        end
        xlabel('time (s)')
        ylabel('Spks/Sec')
        set(gca,'YTick','')
        line([0 0],get(gca,'ylim'),'color',[0 0 0],'LineStyle','--')

        legend([h.mainLine],CondNames)


        % Rayleigh Plots
        figure('name',[thismouse ' ' thisdate ' ' thisses ': Frequency Modulation'],'units','normalized','outerposition',[0 0 1 1])
        for freqid=1:length(freqs)
            cid = find(cell2mat(cellfun(@(X) any(strfind(X,num2str(freqs(freqid)*10))),CondNames,'UniformOutput',0)));
            trialidx = find(ismember(condindx,cid));
            circledur = 1./freqs(freqid);
            angles = arrayfun(@(X) (newtimevec(tmpspks(X,:)>0)*2*pi)./circledur,trialidx,'UniformOutput',0);
            if isempty([angles{:}])
                continue
            end

            subplot(length(freqs),2,(2*freqid)-1)
            histogram([angles{:}]',[0:0.2*pi:max([angles{:}])])
            ylims = get(gca,'ylim');
            arrayfun(@(X) patch([X X+pi X+pi X],[min(ylims) min(ylims) max(ylims) max(ylims)],[0 1 0],'FaceAlpha',0.2,'EdgeColor','none'),0:2*pi:max([angles{:}]),'UniformOutput',0)
            set(gca,'XTickLabel',cellfun(@(X) num2str(round(str2num(X)./(2*pi).*circledur*100)./100),(get(gca,'XTickLabel')),'UniformOutput',0))
            xlabel('Time (s)')
            hold on
            ylabel(['count'])
            title([num2str(freqs(freqid)) 'Hz'])
            makepretty

            subplot(length(freqs),2,2*freqid)
            circ_plot([angles{:}]','hist',[],20,true,true,'linewidth',2,'color','r')
            p=circ_rtest([angles{:}]');
            title(['p=' num2str(p)])
            makepretty
        end



    end
end
