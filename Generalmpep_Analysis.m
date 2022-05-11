%% User Input
% Load all data
% Find available datasets (always using dates as folders)

RedoAfterClustering=1;
Redo =0; % Redo in general
RedoTable=1;
NewHistologyNeeded = 0; %Automatically to 1 after RedoAfterClustering
%Predefine
SaveRFDir = 'E:\Data\Mpep\'
abortsession = 0;
timeBinSize = 1/100;
pretrialtime = 0.2; %take up to x seconds prior trial
posttrialtime = 0.2; % take up to x seconds post trial
% Build wavelet filter bank
freqlims = [0.85 9.79];
nv = 48;
plotthis = 0; %For intermediate step plots (only for debugging)
%% Automated
% Build filterbank for wavelet transforms


clear DateOpt
for idx = 1:length(DataDir)
    DateOpt{idx} = cellfun(@(X) dir(fullfile(DataDir{idx},X,'*-*')),MiceOpt(DataDir2Use==idx),'UniformOutput',0);
end
DateOpt = cat(2,DateOpt{:});
DateOpt = cellfun(@(X) X([X.isdir]),DateOpt,'UniformOutput',0);
DateOpt = cellfun(@(X) {X.name},DateOpt,'UniformOutput',0);
for midx = 1:length(MiceOpt)
    Dates4Mouse = DateOpt{midx};
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
            thisses = mpepsess{sesidx};

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
            
            %% Load Protocol and identify unique conditions
            Protocol = load(fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},thisdate,thisses,'Protocol.mat'));
            Protocol = Protocol.Protocol;
            
            % extract trial information
            ntrials = prod(size(Protocol.seqnums));
            ParFullnames = Protocol.pardefs;
            ParDefs = Protocol.parnames;
            Pars = Protocol.pars;
            ncond = size(Pars,2);
            sequence = Protocol.seqnums;
            [condindx,repidx] = arrayfun(@(X) find(sequence==X),1:ntrials,'UniformOutput',0);
            
            % Factors of relevance
            parrel = [];
            RelName = {};
            for parid = 1:size(Pars,1)
                if length(unique(Pars(parid,:)))>1
                    parrel = [parrel parid];
                    RelName = {RelName{:} ParDefs{parid}};
                end
            end
            
            condindx = cell2mat(condindx); %Condition index in trial ordedr
            TrialDurations = Pars(strcmp(ParDefs,'dur'),condindx)/10;
            
            newtimevec = -pretrialtime:timeBinSize:max(TrialDurations)+posttrialtime;
            timeEdges = -pretrialtime-timeBinSize/2:timeBinSize:max(TrialDurations)+posttrialtime+timeBinSize/2;
            
            %% Align to Trial Onset times
            [starttrialidx,endtrialidx,trialid] = FindTrialOnAndOffSetsPhotoDiode(Actualtime,Timeline(:,ismember(AllInputs,'photoDiode')),TrialDurations);
            
            TrialDurations(isnan(starttrialidx))=[];
            endtrialidx(isnan(starttrialidx))=[];
            starttrialidx(isnan(starttrialidx))=[];
            
            if any((Actualtime(endtrialidx)-Actualtime(starttrialidx)) - TrialDurations'>1/tmSR*10)
                disp(['flips slightly drifting... Average of ' num2str(nanmean((Actualtime(endtrialidx)-Actualtime(starttrialidx)) - TrialDurations')) 'sec'])
                
            end
            if ntrials ~= length(starttrialidx)
                warning('Can''t find enough trials')
                keyboard
            end
            
            %% Loading data from kilosort/phy easily
            myKsDir = fullfile(LocalDir,MiceOpt{midx});
            subksdirs = dir(fullfile(myKsDir,'Probe*')); %This changed because now I suddenly had 2 probes per recording
            
            if strcmp(ProbeType{midx},'2_4S') && ~isempty(subksdirs) % These are my chronic mice, one dataset per mouse
              myKsDir = fullfile(LocalDir,MiceOpt{midx});
                subksdirs = dir(fullfile(myKsDir,'Probe*')); %This changed because now I suddenly had 2 probes per recording
                if length(subksdirs)<1
                    clear subksdirs
                    subksdirs.folder = myKsDir; %Should be a struct array
                    subksdirs.name = 'Probe0';
                end
            else
                myKsDir = fullfile(LocalDir,MiceOpt{midx},thisdate);
                subksdirs = dir(fullfile(myKsDir,'Probe*')); %This changed because now I suddenly had 2 probes per recording
                if length(subksdirs)<1
                    clear subksdirs
                    subksdirs.folder = myKsDir; %Should be a struct array
                    subksdirs.name = 'Probe0';
                end
            end
            
            for probeid = 1:length(subksdirs)
               
                myKsDir = fullfile(subksdirs(probeid).folder,subksdirs(probeid).name);
                if ~isdir(myKsDir)
                    continue
                end
                
                %Saving directory
                thisprobe = subksdirs(probeid).name;
                if ~Redo && exist(fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,'MPEPData.mat'))
                    disp([MiceOpt{midx} ' ' thisdate   ' ' thisses ' already processed... skipping'])
                    continue
                end
                
                if ~Redo && exist(fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,'RFData.mat'))
                    if ~RedoAfterClustering || exist(fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,thisprobe,'CuratedResults.mat'))
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
                
                %Saving directory
                if ~isdir(fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe))
                    mkdir(fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe))
                end
                
                %% Computing some useful details about spikes/neurons (like depths)
                myLFDir = fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},thisdate,'ephys');
                lfpD = dir(fullfile(myLFDir,'*','*','*.ap.*bin')); % ap file from spikeGLX specifically
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
                clear params
                params.loadPCs=true;
                params.thisdate = thisdate;
                PrepareClusInfo
                
                %% Get Histology output
                if strcmp(ProbeType{midx},'2_4S')
                    thisdate = []; % There's no data for the chronic mice in front of histology.
                end
                GetHistologyOutput      %I created an extra function to have one line of code in the different scripts
                if ~histoflag
                    disp([MiceOpt{midx} thisdate thisses thisprobe 'No histology data...'])
                end
                thisdate = Dates4Mouse{didx}; % Reassign thisdate
                
                %% load synchronization data
                SyncKSDataToTimeline
                
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
                nclus = length(Good_IDx);
                if nclus < 2
                    disp('Less than 2 good units, skip...')
                    continue
                end
                
                %% Spikes per trial
                SpikeIDx = arrayfun(@(X) find(spikeTimesCorrected>=Actualtime(starttrialidx(X))-pretrialtime&spikeTimesCorrected<=Actualtime(endtrialidx(X))+posttrialtime),1:ntrials,'UniformOutput',0);
                SpikeTrialID = nan(1,length(spikeTimesCorrected));
                SpikeTrialTime = nan(1,length(spikeTimesCorrected));
                for trid = 1:ntrials
                    SpikeTrialID(SpikeIDx{trid})=trid;
                    SpikeTrialTime(SpikeIDx{trid}) = spikeTimesCorrected(SpikeIDx{trid})-Actualtime(starttrialidx(trid));
                end
                %remove NaNs
                SpikeTrialTime(isnan(SpikeTrialID))=[];
                spikeCluster(isnan(SpikeTrialID)) = [];
                spikeTimesCorrected(isnan(SpikeTrialID))=[];
                SpikeTrialID(isnan(SpikeTrialID))=[];
                
                %Spike rate / histogram
                SpikeRatePerTP = arrayfun(@(Y) arrayfun(@(X) histcounts(SpikeTrialTime(SpikeTrialID == X & spikeCluster'== Y),...
                    timeEdges),1:ntrials,'UniformOutput',0),cluster_id(Good_IDx),'UniformOutput',0);
                SpikeRatePerTP = cat(1,SpikeRatePerTP{:});
                SpikeRatePerTP = cat(1,SpikeRatePerTP{:});
                SpikeRatePerTP = reshape(SpikeRatePerTP,length(Good_IDx),ntrials,[]); %Reshape to nclus, ntrials, ntp
                SpikeRatePerTP=SpikeRatePerTP./timeBinSize; %in spikes/sec
                SpikeRatePerTP = permute(SpikeRatePerTP,[2,3,1]); % Convert to trial,time,unit
                SpikeRatePerTP(repmat(sum(SpikeRatePerTP==0,3)==nclus,[1,1,nclus]))=nan;% Fill with nans instead of 0 when longer
                if ~any(SpikeRatePerTP(:)>0)
                    disp(['No spikes for ' MiceOpt{midx} ' ' thisdate ' ' thisses ', skip..'])
                    continue
                end
                %% Plot per condition
                figure('name','Average PSTH')
                [condsorted,sortid] = sort(condindx);
                imagesc(newtimevec,condsorted,nanmean(SpikeRatePerTP(sortid,:,:),3))
                colormap hot
                ylabel('Condition')
                xlabel('Time (s)')
                saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,['AveragePSTH.fig']))
                saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,['AveragePSTH.bmp']))
                
                cols = jet(ncond);
                
                figure('name','PSTH per condition')
                AllCondNames = {};
                for condid = 1:ncond
                    subplot(ceil(sqrt(ncond)),round(sqrt(ncond)),condid)
                    tmp = squeeze(nanmean(SpikeRatePerTP(condindx==condid,:,:),1));
                    h=shadedErrorBar(newtimevec,nanmean(tmp,2),nanstd(tmp,[],2)./sqrt(nclus-1),'lineProps',{'color',cols(condid,:),'LineWidth',1.5});
                    ylabel('Spikes/sec')
                    xlabel('Time (s)')
                    conditionvals = arrayfun(@(X) [' ' ParDefs{X} '=' num2str(Pars(X,condid))],parrel,'UniformOutput',0);
                    conditionvals = {[conditionvals{:}]};
                    AllCondNames = {AllCondNames{:} conditionvals{1}};
                    title(conditionvals)
                end
                saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,['AveragePSTHPerCondition.fig']))
                saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,['AveragePSTHPerCondition.bmp']))
                
                %% TF analysis
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
                figure('name','Average amplitude across trials');
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
                
                figure('name','TF per condition')
                for condid = 1:ncond
                    subplot(ceil(sqrt(ncond)),round(sqrt(ncond)),condid)
                    % Spikes
                    tmpspks = nanmean(nanmean(SpikeRatePerTP(condindx==condid,:,:),1),3);
                    tmpspks(isnan(tmpspks))=0;
                    [csf,f] = cwt(double(tmpspks),'FilterBank',fb);
                    Ampl= abs(csf);
                    
                    
                    h=imagesc(newtimevec,[],Ampl);
                    xlabel('Time (s)')
                    ylabel('Frequency Indx')
                    colormap hot
                    title(AllCondNames{condid})
                    hold on
                    
                    yyaxis right
                    avgAmpl = nanmean(Ampl,2);
                    plot(avgAmpl,f,'b-');
                    ylabel('Frequency (Hz)')
                    set(gca,'YScale','log')
                    makepretty
                end
                saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,['TFperCondition.fig']))
                saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,['TFperCondition.bmp']))
                
                
                
                %% Initialize save data
                MpepInfo = table;
                MpepInfo.ClusID = cluster_id(Good_IDx);
                MpepInfo.depth = depth(Good_IDx);
                MpepInfo.SpikesPerSec = permute(SpikeRatePerTP,[3,1,2]);
                
                
                [sorteddepth,sortidx] = sort(depth(Good_IDx));
                
                TrialInfo = table;
                TrialInfo.CondNames = arrayfun(@(X) AllCondNames{X},condindx,'UniformOutput',0)';
                TrialInfo.CondIndx = condindx';
                %% Analyze without plotting for all units
                if any(ismember(RelName,{'tf','ModFreq'})) % Frequency modulation!
                    parvec = find(ismember(ParDefs,{'tf','ModFreq'}));
                    
                    figure('name',['Freq modulation across depth - Rayleigh']);
                    
                    for parid = 1:length(parvec)
                        freqs = unique(Pars(parvec(parid),:))/10;
                        
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
                        eval(['MpepInfo.' ParDefs{parvec(parid)} '_Rvec=Rvec;'])
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
                        areatmp = Depth2AreaPerUnit.Area(Good_IDx);
                        areatmp = strrep(areatmp,'/','');
                        areacol = cellfun(@(X) hex2rgb(X),Depth2AreaPerUnit.Color(Good_IDx),'UniformOutput',0);
                        
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
                
                %% Plot per unit
                nexample = 5
                
                %                 exampleunits = datasample(find(any(RayleighP<0.001,2)),nexample,'replace',false);
                exampleunits = datasample(1:nclus,nexample,'replace',false);
                
                for uid = 1:nexample
                    unitid = cluster_id(Good_IDx(exampleunits(uid)));
                    depthhere = depth(Good_IDx(exampleunits(uid)));
                    if histoflag
                        unitname = ['Unit ' num2str(unitid) ', depth=' num2str(depthhere) ', ' areatmp{exampleunits(uid)}];
                    else
                        unitname = ['Unit ' num2str(unitid) ', depth=' num2str(depthhere)];
                    end
                    
                    figure('name',unitname)
                    subplot(2,1,1)
                    hold on
                    arrayfun(@(X) scatter(newtimevec(SpikeRatePerTP(sortid(X),:,exampleunits(uid))>0),repmat(X,[1,sum(SpikeRatePerTP(sortid(X),:,exampleunits(uid))>0)]),8,cols(condsorted(X),:),'filled'),1:ntrials,'UniformOutput',0)
                    xlabel('Time (s)')
                    ylabel('Trial (sorted by condition)')
                    line([0 0],get(gca,'ylim'),'color',[0 0 0],'LineStyle','--')
                    
                    subplot(2,1,2)
                    hold on
                    clear h
                    for condid=1:ncond
                        tmp = SpikeRatePerTP(condindx==condid,:,exampleunits(uid));
                        h(condid)=shadedErrorBar(newtimevec,smooth(nanmean(tmp,1)+max(get(gca,'ylim'))*0.9),smooth(nanstd(tmp,[],1)./sqrt(size(tmp,1)-1)),'lineProps',{'color',cols(condid,:),'LineWidth',1.5});
                    end
                    xlabel('time (s)')
                    ylabel('Spks/Sec')
                    set(gca,'YTick','')
                    line([0 0],get(gca,'ylim'),'color',[0 0 0],'LineStyle','--')
                    
                    legend([h.mainLine],AllCondNames)
                    saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,[unitname '_PSTH.fig']))
                    saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,[unitname '_PSTH.bmp']))
                    
                    
                    if any(ismember(RelName,{'tf','ModFreq'})) % Frequency modulation!
                        parvec = find(ismember(ParDefs,{'tf','ModFreq'}));
                        for parid = 1:length(parvec)
                            freqs = unique(Pars(parvec(parid),:))/10;
                            
                            
                            figure('name',[unitname ParDefs{parvec(parid)} ' Frequency Modulation'])
                            for freqid=1:length(freqs)
                                
                                trialidx = find(ismember(condindx,find(Pars(parvec(parid),:)==freqs(freqid)*10)));
                                circledur = 1./freqs(freqid);
                                angles = arrayfun(@(X) (newtimevec(SpikeRatePerTP(X,:,exampleunits(uid))>0)*2*pi)./circledur,trialidx,'UniformOutput',0);
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
                                ylabel(['nr spikes at phase ' ParDefs{parvec(parid)}])
                                title([num2str(freqs(freqid)) 'Hz'])
                                makepretty
                                
                                subplot(length(freqs),2,2*freqid)
                                circ_plot([angles{:}]','hist',[],20,true,true,'linewidth',2,'color','r')
                                p=circ_rtest([angles{:}]');
                                title(['p=' num2str(p)])
                                makepretty
                            end
                            
                            saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,[unitname ParDefs{parvec(parid)} ' Frequency Modulation.fig']))
                            saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,[unitname ParDefs{parvec(parid)} ' Frequency Modulation.bmp']))
                            
                        end
                        
                    end
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


%% table across mice
if ~RedoTable || ~exist(fullfile(SaveDir,'AllMiceMPEPData.mat'))
    AllMPEPDat = [];
    AllTrialInfo = table;
    for midx = 1:length(MiceOpt)
        Dates4Mouse = DateOpt{midx};
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
                    TrialInfo.Probe = repmat(['Probe' ProbeId{2}],ntrials,1);
                   
                    % Add area to mpepinfo table
                    if isfield(tmp,'Depth2AreaPerUnit')
                    tmp.MpepInfo.Area = tmp.Depth2AreaPerUnit.Area(ismember(tmp.Depth2AreaPerUnit.Cluster_ID,tmp.MpepInfo.ClusID));
                    tmp.MpepInfo.Area  = strrep(tmp.MpepInfo.Area,'/','');
                    
                    tmp.MpepInfo.Color = tmp.Depth2AreaPerUnit.Color(ismember(tmp.Depth2AreaPerUnit.Cluster_ID,tmp.MpepInfo.ClusID));
                    end
                    % Add session&mouse&probe information
                    tmp.MpepInfo.Mouse = repmat(MiceOpt{midx},length(tmp.MpepInfo.ClusID),1);
                    tmp.MpepInfo.Date = repmat(thisdate,length(tmp.MpepInfo.ClusID),1);
                    tmp.MpepInfo.Session = repmat(thisses,length(tmp.MpepInfo.ClusID),1);
                    tmp.MpepInfo.Probe = repmat(thisprobe,length(tmp.MpepInfo.ClusID),1);

                    
                    % Save in big struct
                    if isempty(AllMPEPDat)                        
                        AllMPEPDat = tmp.MpepInfo;
                        AllMPEPDat = tall(AllMPEPDat); %Make tall to prevent memory issues                        
                        AllTrialInfo = TrialInfo;
                    else
                        % Concatenate
                        AllMPEPDat = cat(1,AllMPEPDat,tmp.MpepInfo);
                        AllTrialInfo = cat(1,AllTrialInfo,TrialInfo);                        
                    end
                end
            end
        end
    end
    save(fullfile(SaveRFDir,'AllMiceMPEPData.mat'),'AllMPEPDat','AllTrialInfo','-v7.3')
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
figure('name','VisualVsAuditorymodulation')
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
figure('name','Pvalue Distribution')
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
figure('name','Significant Modulation')
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
    set(gca,'XTickLabel',AutomaticAREASOfInterest(include),'XTickLabelRotation',90)
    ylabel('Nr Units')
    linkaxes([ax1 ax2])
end

%% Example scatter
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
condid = AllTrialInfo.CondIndx;
CondNames = AllTrialInfo.CondNames;
[condopt,idx,udx] = unique(condid);
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
    
    figure('name',[thismouse ' ' thisdate ' ' thisses ': Unit ' num2str(clusterid(exampleunits(uid))) areatmp{exampleunits(uid)}])

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
    for condid=1:ncond
        tmp = tmpspks(condindx==condid,:);
        h(condid)=shadedErrorBar(newtimevec,smooth(nanmean(tmp,1)+max(get(gca,'ylim'))*0.9),smooth(nanstd(tmp,[],1)./sqrt(size(tmp,1)-1)),'lineProps',{'color',cols(condid,:),'LineWidth',1.5});
    end
    xlabel('time (s)')
    ylabel('Spks/Sec')
    set(gca,'YTick','')
    line([0 0],get(gca,'ylim'),'color',[0 0 0],'LineStyle','--')
    
    legend([h.mainLine],CondNames)
    
    
    % Rayleigh Plots
    figure('name',[thismouse ' ' thisdate ' ' thisses ': Frequency Modulation'])
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
