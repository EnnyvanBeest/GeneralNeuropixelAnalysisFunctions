%% User Input
% Load all data
% Find available datasets (always using dates as folders)
usehistology = 0; %if 0, just use depths, can actually help histology alignment later
Redo=0; % Redo in general
RedoAfterClustering = 0; %keep 0,
%Predefine
SaveRFDir = SaveDir
abortsession = 0;
nboot = 100;
drawthis = 1;
clear Depth2AreaPerUnit
%% Automated
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
        RFsess = [];
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
                
                if any(strfind(A,'SparseNoise'))
                    % Load protocol to check for contrast >0
                    try
                        Protocol = load(fullfile(fullfile(listfiles(idx).folder,'Protocol.mat')));
                    catch ME
                        disp(ME)
                        continue
                    end
                    Protocol = Protocol.Protocol;
                    if Protocol.pars(strcmp(Protocol.parnames,'c')) == 0
                        %not a RF session, continue
                        continue
                    end
                    
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
        
        for sesidx = 1:length(RFsess)
            close all
            abortsession = 0;
            thisses = RFsess{sesidx};
            
            %% which probes?
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
            
            % Check if exists already, and if so redo?
            Runthis = 0;
            for probeid = 1:length(subksdirs)
                histoflag = 0;
                NewHistologyNeeded = 0;
                myKsDir = fullfile(subksdirs(probeid).folder,subksdirs(probeid).name);
                if ~isdir(myKsDir)
                    continue
                end
                
                %Saving directory
                thisprobe = subksdirs(probeid).name;
                if ~Redo && exist(fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,'RFData.mat'))
                    if ~RedoAfterClustering || exist(fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,thisprobe,'CuratedResults.mat'))
                        continue
                    elseif RedoAfterClustering
                        myKsDir = fullfile(LocalDir,MiceOpt{midx},thisdate,thisprobe);
                        myClusFile = dir(fullfile(myKsDir,'cluster_info.tsv'));
                        if isempty(myClusFile)
                            disp('This data is not yet curated with phy!!')
                            continue
                        end
                        NewHistologyNeeded = 1; %Automatically to 1 after RedoAfterClustering
                    end
                else
                    Runthis=1;
                end
            end
            if ~Runthis
                disp(['No need to waste time on this session ' MiceOpt{midx} ' ' thisdate ' ' thisses])
                continue
            end
            %% Timeline
            disp(['Running ' MiceOpt{midx} ' ' thisdate ' ' thisses])
            
            % Timeline to go to:
            timelinefile = dir(fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},thisdate,thisses,'*Timeline.mat'));
            timelineabsentwarning =0;
            if isempty(timelinefile)
                timelineabsentwarning =0;
                disp([MiceOpt{midx} ' ' thisdate   ' ' thisses ' has no timeline ... skipping'])
                continue
            end
            loadTimeline
            
            %% Load RF mapping session details
            if exist(fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},thisdate,thisses,'RF_MappingStimuli.mat'))
                SS = load(fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},thisdate,thisses,'RF_MappingStimuli.mat'));
            else
                disp([fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},thisdate,thisses,'RF_MappingStimuli.mat') ' does not yet exist.. create on rig using RFmapping_RegenerateStimuli'])
                continue
            end
            if isempty(SS)
                continue
            end
            
            Protocol = load(fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},thisdate,thisses,'Protocol.mat'));
            Protocol = Protocol.Protocol;
            % extract trial information
            SS = SS.SS;
            ntrials = length(SS);
            TrialDurations = cell2mat(cellfun(@(X) length(X.ImageSequence)./60,SS,'UniformOutput',0)); %in sec*10
            visSpace = [0 1; 0 1]; % Proportion of azimuth, proportion of elevation, 0 = left/top, 1 = right/bottom of stimulus.
            
            %% Align to Trial Onset times
            [starttrialidx,endtrialidx,trialid] = FindTrialOnAndOffSetsPhotoDiode(Actualtime,Timeline(:,ismember(AllInputs,'photoDiode')),TrialDurations);
            
            TrialDurations(isnan(starttrialidx))=[];
            endtrialidx(isnan(starttrialidx))=[];
            starttrialidx(isnan(starttrialidx))=[];
            ntrials = length(starttrialidx);
            if ntrials==0
               disp('Cant find proper data alignment.. skip')
                continue
            end
            if any((Actualtime(endtrialidx)-Actualtime(starttrialidx)) - TrialDurations'>1/tmSR*10)
                disp(['flips slightly drifting... Average of ' num2str(nanmean((Actualtime(endtrialidx)-Actualtime(starttrialidx)) - TrialDurations')) 'sec'])
            end
            if ntrials ~= length(starttrialidx)
                warning('Can''t find enough trials')
                keyboard
            end
            
            %% Prepare Noise stimulus frames in time X rows X Cols
            tmp = Timeline(:,ismember(AllInputs,'photoDiode'));
            %             tmp = detrend(tmp);
            tmp = (tmp-nanmin(tmp(:)))./(nanmax(tmp(:))-nanmin(tmp(:)));
            tresh=nanmedian(tmp);
            tmp(tmp>=tresh)=1;
            tmp(tmp<tresh) = 0;
            photodiodetrace = tmp;
            
            figure; plot(Actualtime,photodiodetrace);
            box off
            tmpimg = cat(3,SS{1}.ImageTextures{:});
            uniquecolors = unique(tmpimg(:));
            Greyvalue = uniquecolors(2);
            stimFrames = repmat(Greyvalue,ceil(length(Actualtime)./(tmSR/60)),size(tmpimg,1),size(tmpimg,2));
            countimagetotal = 1;
            Stimulustimes = nan(1,ceil(length(Actualtime)./(tmSR/60)));
            for tridx=1:ntrials
                countimg = 0;
                counttime = 1/60;
                diodestatus = photodiodetrace(1);
                missedflips = [];
                framesmissing = 0;
                for tp=starttrialidx(tridx):endtrialidx(tridx)
                    % photodiodestate
                    if photodiodetrace(tp)~= diodestatus %60 Hz flips
                        %                         if counttime>1/60*1.6
                        %                             xlim([Actualtime(tp-500) Actualtime(tp+500)])
                        %                             if exist('h')
                        %                                 delete(h)
                        %                             end
                        %                             h=line([Actualtime(tp) Actualtime(tp)],get(gca,'ylim'),'color',[1 0 0],'LineStyle','--');
                        %                             disp('too slow')
                        %                             missedflips = [missedflips tp];
                        %                             drawnow
                        %                             pause(0.5)
                        %
                        %                         end
                        if counttime<1/60*0.6
                            xlim([Actualtime(tp-500) Actualtime(tp+500)])
                            if exist('h')
                                delete(h)
                            end
                            h=line([Actualtime(tp) Actualtime(tp)],get(gca,'ylim'),'color',[1 0 0],'LineStyle','--');
                            drawnow
                            %                             pause(0.5)
                        else
                            countimg = countimg+1; %go to next image
                            %Load current image
                            if countimg > length(SS{tridx}.ImageTextures)
                                xlim([Actualtime(tp-500) Actualtime(tp+500)])
                                if exist('h')
                                    delete(h)
                                end
                                h=line([Actualtime(tp) Actualtime(tp)],get(gca,'ylim'),'color',[1 0 0],'LineStyle','--');
                                
                                
                                countimg = length(SS{tridx}.ImageTextures);
                                framesmissing = framesmissing+1;
                                drawnow
                                pause(0.5)
                                
                                if length(framesmissing)>1 %1 last flip is grey screen (should be)
                                    disp('More Photodiode flips than frames?')
                                    keyboard
                                end
                            end
                            
                            tmpimg =  SS{tridx}.ImageTextures{SS{tridx}.ImageSequence(countimg)};
                            stimFrames(countimagetotal,:,:) = tmpimg;
                            Stimulustimes(countimagetotal)=Actualtime(tp);
                            countimagetotal = countimagetotal+1;
                            
                            counttime = 0;
                            diodestatus=photodiodetrace(tp);
                        end
                    end
                    counttime = counttime+(1/tmSR);
                end
            end
            %interpolate nans stimulustime
            Stimulustimes = fillmissing(Stimulustimes,'linear');
            disp(['Individual square are shown for '  num2str(length(Stimulustimes)./countimagetotal)])
            
            nX = size(stimFrames,3);
            nY = size(stimFrames,2);
            xPos = 1:nX;
            yPos = 1:nY;
            %                 crosscorrelation to see how much stimuli
            %                 overlap
            xcor= corr(reshape(stimFrames,[],nY*nX));
            figure; imagesc(xcor,[-0.5 0.5]); title('Sanity check: Correlation in stimulus presentations')
            
            stimFramesOri=stimFrames; %When having >1 probes this is important
            
            stimFrameDur = median(diff(Stimulustimes));
            tmptime = downsample(Actualtime,15);
            newsr = nanmedian(diff(tmptime));
            newtimevec = Stimulustimes(1):newsr:Stimulustimes(end);
            newtimebins = [newtimevec(1)-newsr/2:newsr:newtimevec(end)+newsr/2];
            
            %% Loading data from kilosort/phy easily
            for probeid = 1:length(subksdirs)
                histoflag = 0;
                NewHistologyNeeded = 0;
                myKsDir = fullfile(subksdirs(probeid).folder,subksdirs(probeid).name);
                if ~isdir(myKsDir)
                    continue
                end
                
                %Saving directory
                thisprobe = subksdirs(probeid).name;
                if ~Redo && exist(fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,'RFData.mat'))
                    if ~RedoAfterClustering || exist(fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,thisprobe,'CuratedResults.mat'))
                        continue
                    elseif RedoAfterClustering
                        myKsDir = fullfile(LocalDir,MiceOpt{midx},thisdate,thisprobe);
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
                else
                    % move everything in this directory to subdirectory
                    % with 'old'
%                     movefile(fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe),fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,'old'))
                end
                
                %% Computing some useful details about spikes/neurons (like depths)
                myLFDir = fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},thisdate,'ephys');
                lfpD = dir(fullfile(myLFDir,'*','*','*.ap.*bin')); % ap file from spikeGLX specifically
                if length(lfpD)<length(subksdirs)
                    lfpD = dir(fullfile(myLFDir,'*VR*','*','*.ap.*bin')); % ap file from spikeGLX specifically
                    if length(lfpD)~=length(subksdirs)
                        disp('Should be a different amount of probes?')
                        keyboard
                    end
                end
                lfpD = lfpD(probeid);
                
                % Get information from meta file
                [Imecmeta] = ReadMeta2(lfpD.folder,'ap');
                lfpFs = str2num(Imecmeta.imSampRate);
                imecSR = 1./lfpFs;
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
                if usehistology
                    GetHistologyOutput      %I created an extra function to have one line of code in the different scripts
                    if ~histoflag
                        disp([MiceOpt{midx} thisdate thisses thisprobe 'No histology data...'])
                    end
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

                %% Divide spike times in trials and clusters

                SpikeCountPerTP = arrayfun(@(X) histcounts(spikeTimesCorrected(spikeCluster'==X),newtimebins),cluster_id(Good_IDx),'UniformOutput',0);
                SpikeCountPerTP = cat(1,SpikeCountPerTP{:})';
%                 % Prepare spikecounts per unit according to actualtime X units
%                 SpikeCountPerTP = zeros(length(newtimevec),length(Good_IDx));
%                 parfor clusid = 1:length(Good_IDx)
%                     clusid
%                     SpikesThisCluster = spikeTimesCorrected(spikeCluster'== cluster_id(Good_IDx(clusid)));
%                     if isempty(SpikesThisCluster)
%                         continue
%                     end
%                     % spikecount per stimulus presentation time point,
%                     SpikeCountPerTP(:,clusid) = histcounts(SpikesThisCluster,newtimebins);
%                 end
                
                %RunningSpeed
                RotarySignal = (Timeline(:,ismember(AllInputs,'rotaryEncoder')));
                RotarySignalPerSec = nan(1,length(newtimevec));
                windowsz = tmSR; %per second
                
                parfor tp = 1:length(newtimevec)
                    RotarySignalPerSec(tp) = nanmean(diff(RotarySignal(Actualtime>=newtimebins(tp)&Actualtime<=newtimebins(tp+1))));
                end
                % Threshold
                RotarySignalPerSec(abs(RotarySignalPerSec)>nanmedian(abs(RotarySignalPerSec))) = 1;
                %% Make space
                clear tmp
                clear tmp2
                clear tmpdat
                clear tmpdiff
                clear tmpspikeshere
                clear SS
                clear RotarySignal
                clear photodiode
                clear N
                clear photodiodetrace
                
                
                %% Prepare data for fit
                theta = linspace(0,2*pi,100);
                
                options = optimset('fminsearch');
                options.MaxFunEvals = 1000;
                options.Display = 'iter';
                
                RFSizeReq = 5; %Required in visual degree for analysis
                RFSizeCur = Protocol.pars(7)/10;
                %                     Resample images by that much
                RSize = RFSizeCur./RFSizeReq;
                if RSize>1
                    RSize=1;
                    RFSizeReq = RFSizeCur;
                end
                
                stimFrames = permute(stimFramesOri,[2,3,1]); % Shift dimensions for imresize
                stimFrames = imresize(stimFrames,RSize); % Resize image
                newY = size(stimFrames,1); %New dimensions
                newX = size(stimFrames,2);
                stimFrames = permute(stimFrames,[3,1,2]); % Dimensions back to ori
                % For every receptive field make a trace with 1 for
                % white and -1 for black pixel
                stimFrames = reshape(stimFrames,[],newY*newX);
                greyval = nanmedian(unique(stimFrames(:)));
                whiteval = quantile(unique(stimFrames(:)),0.75);
                blackval = quantile(unique(stimFrames(:)),0.25);
                stimFrames(stimFrames<blackval)=-1;
                stimFrames(stimFrames>whiteval)=1;
                stimFrames(abs(stimFrames)~=1)=0;
                
                % put in same timeline
                NewStimVecON = zeros(length(newtimevec),newY*newX);
                NewStimVecOFF = zeros(length(newtimevec),newY*newX);
                parfor tp = 1:length(newtimevec)
                    tptotake = Stimulustimes>=newtimebins(tp)&Stimulustimes<=newtimebins(tp+1);
                    if any(tptotake)
                        NewStimVecON(tp,:) = sign(nanmean(stimFrames(tptotake,:),1))==1;
                        NewStimVecOFF(tp,:) = sign(nanmean(stimFrames(tptotake,:),1))==-1;
                    end
                end
                
                %                 InclTP = find(sum(SpikeCountPerTP==0,2)<0.8*size(SpikeCountPerTP,2) & ~isnan(RotarySignalPerSec)');
                InclTP = 1:size(SpikeCountPerTP,1);%find(sum(SpikeCountPerTP==0,2)<0.8*size(SpikeCountPerTP,2) & ~isnan(RotarySignalPerSec)');
                InclUnits = find(sum(SpikeCountPerTP,1)>10); %need at least 10 spikes
                
                if isempty(InclUnits)
                    disp(['No Units left for ' MiceOpt{midx} ', ' thisdate ', ' thisses])
                    continue
                end
                okay=0;
                breaksession = 0;
                while ~okay
                    % Regress out running
                    %                     DM = cat(2,ones(length(InclTP),1),RotarySignalPerSec(InclTP)');
                    
                    %GLM
                    if histoflag
                        areashere = Depth2AreaPerUnit.Area(ismember(Depth2AreaPerUnit.Cluster_ID,cluster_id(Good_IDx(InclUnits))));
                        areashere = cellfun(@(X) strrep(X,'/',''),areashere,'UniformOutput',0)
                    else
                        areashere = arrayfun(@(X) [num2str(depth(X)) 'um'],Good_IDx(InclUnits),'UniformOutput',0);
                    end
                    
                    
                    % Convert to visual degrees
                    xwidth = linspace(Protocol.pars(2)/10,Protocol.pars(3)/10,newX);
                    yheight = linspace(Protocol.pars(4)/10,Protocol.pars(5)/10,newY);
                    
                    
                    %% Evoked RFs
                    
                    % EVoked visual responses
                    tmp = SpikeCountPerTP(InclTP,InclUnits);
                    
                    RFON = arrayfun(@(X) nanmean(tmp(NewStimVecON(InclTP,X)==1,:),1),1:newY*newX,'UniformOutput',0);
                    RFON = cat(1,RFON{:});
                    RFOFF = arrayfun(@(X) nanmean(tmp(NewStimVecOFF(InclTP,X)==1,:),1),1:newY*newX,'UniformOutput',0);
                    RFOFF = cat(1,RFOFF{:});
                    
                    % Again, Z-score
                    RFON = (RFON - nanmean(RFON,1))./nanstd(RFON,[],1);
                    RFOFF = (RFOFF - nanmean(RFOFF,1))./nanstd(RFOFF,[],1);
                    
                    
                    %Smooth
                    %                     RFAll = arrayfun(@(X) smooth2a(reshape(RFAll(:,X),newY,newX),2),1:size(RFAll,2),'UniformOutput',0);
                    %                     RFAll = cat(3,RFAll{:});
                    RFON = arrayfun(@(X) smooth2a(reshape(RFON(:,X),newY,newX),2),1:size(RFON,2),'UniformOutput',0);
                    RFON = cat(3,RFON{:});
                    RFOFF = arrayfun(@(X) smooth2a(reshape(RFOFF(:,X),newY,newX),2),1:size(RFOFF,2),'UniformOutput',0);
                    RFOFF = cat(3,RFOFF{:});
                    
                    RFAll = abs(RFON) + abs(RFOFF);
                    
                    if sum(sum(isnan(nanmean(RFAll,3))))>0.25*newY*newX
                        disp('Not enough data.. skip session')
                        breaksession = 1;
                        
                        break
                    end
                    % Interpollate nans
                    RFAll = fillmissing(RFAll,'spline');
                    RFON = fillmissing(RFON,'spline');
                    RFOFF = fillmissing(RFOFF,'spline');
                    if any(any(isnan(nanmean(RFAll,3))))
                        disp('Not enough data.. skip session')
                        breaksession = 1;
                        break
                    end
                    if any(any(isnan(reshape(RFAll,newY*newX,[])),1))
                        InclUnits(any(isnan(reshape(RFAll,newY*newX,[])),1)) = [];
                    else
                        okay = 1;
                    end
                end

                %% Z-score
                RFAll = reshape(RFAll,newY*newX,[]);
                RFAll = reshape((RFAll-nanmean(RFAll,1))./nanstd(RFAll,[],1),newY,newX,[]);

                RFON = reshape(RFON,newY*newX,[]);
                RFON = reshape((RFON-nanmean(RFON,1))./nanstd(RFON,[],1),newY,newX,[]);

                RFOFF = reshape(RFOFF,newY*newX,[]);
                RFOFF = reshape((RFOFF-nanmean(RFOFF,1))./nanstd(RFOFF,[],1),newY,newX,[]);
                if breaksession
                    continue
                end
                
                %% Fit 2D GaussianParameters
                options = optimoptions('lsqcurvefit');
                options.Display = 'off';
                options.OptimalityTolerance = 10.^-12;
                options.StepTolerance = 10.^-10;
                [X,Y]=meshgrid(xwidth,yheight);
                lb = [0,min(xwidth),5,min(yheight),5,-pi/4]; %Expect a positive peak!
                ub = [max([RFON(:);RFOFF(:);RFAll(:)]),max(xwidth),max(xwidth)/3,max(yheight),max(yheight)/3,pi/4];
                x0 = [1,0,15,0,15,0]; % Initial guess [Amplitude, x0, sigmax, y0, sigmay, angel(rad)]
                
                %% Fit best parameters to all clusters
                x0 = repmat(x0,[length(InclUnits),1]);
                [x0(:,1),id] = nanmax(reshape(RFAll,[],length(InclUnits)),[],1); % Initial guess [Amplitude, x0, sigmax, y0, sigmay, angel(rad)]
                x0(:,2) = X(id);
                x0(:,4) = Y(id);
                [thisparamsAll,resnormAll,~,~] = arrayfun(@(Z) lsqcurvefit(@D2GaussFunctionRot,x0(Z,:),cat(3,X,Y),RFAll(:,:,Z),lb,ub,options),...
                    1:length(InclUnits),'UniformOutput',0);
                thisparamsAll = cat(1,thisparamsAll{:});
                resnormAll = cat(1,resnormAll{:});
                % predicted values
                       % residual
%                 residual = reshape(RFAll-pred,[],length(InclUnits));%./nanmax(reshape(pred,[],length(InclUnits)),[],1); %Divide by amplitude (to make sure higher amplitudes in the model are not punished!)
                % RFs
%                 resnormAll2 = sum(residual.^2,1);
                
                x0 = thisparamsAll; %replace initial parameters to save time
                [thisparamsON,resnormON,~,~] = arrayfun(@(Z) lsqcurvefit(@D2GaussFunctionRot,x0(Z,:),cat(3,X,Y),RFON(:,:,Z),lb,ub,options),...
                    1:length(InclUnits),'UniformOutput',0);
                thisparamsON = cat(1,thisparamsON{:});
                resnormON = cat(1,resnormON{:});

                % predicted values
                    % residual
%                 residual = reshape(RFON-pred,[],length(InclUnits));%./nanmax(reshape(pred,[],length(InclUnits)),[],1); %Divide by amplitude (to make sure higher amplitudes in the model are not punished!)
%                 resnormON = sum(residual.^2,1);
%                 
                [thisparamsOFF,resnormOFF,~,~] = arrayfun(@(Z) lsqcurvefit(@D2GaussFunctionRot,x0(Z,:),cat(3,X,Y),RFOFF(:,:,Z),lb,ub,options),...
                    1:length(InclUnits),'UniformOutput',0);
                thisparamsOFF = cat(1,thisparamsOFF{:});
                resnormOFF = cat(1,resnormOFF{:});

                % predicted values
%                 resnormOFF = sum(residual.^2,1);
%                 
                
%                 
%                 %                 unitid=95
%                                 thisparams = thisparamsboot;
%                                 figure;
%                                 subplot(2,2,1)
%                                 imagesc(xwidth,yheight,pred(:,:,unitid))
%                                 hold on
%                                 ellipse(thisparams(unitid,3),thisparams(unitid,5),-thisparams(unitid,6),thisparams(unitid,2),thisparams(unitid,4),'g')
%                                 title('predicted')
%                 
%                                 subplot(2,2,2)
%                                 imagesc(xwidth,yheight,RFAll(:,:,unitid))
%                                 hold on
%                                 ellipse(thisparams(unitid,3),thisparams(unitid,5),-thisparams(unitid,6),thisparams(unitid,2),thisparams(unitid,4),'g')
%                                 title('actual')
%                 
%                                 subplot(2,2,3)
%                                 imagesc(xwidth,yheight,reshape(residual(:,unitid).^2,newY,newX))
%                                 hold on
%                                 ellipse(thisparams(unitid,3),thisparams(unitid,5),-thisparams(unitid,6),thisparams(unitid,2),thisparams(unitid,4),'g')
%                                 title('residual')
%                 % %
%                 %
%                 
%                 figure;
%                 for uid=1:length(InclUnits)
%                     imagesc(RFAll(:,:,uid))
%                     pause(0.2)
%                 end
                
                %% Repeat whole thing 1000 times (randomized!) to generate H0-distribution
                resnormbootall = nan(nboot,length(InclUnits)); %all, on, off
                resnormbooton = nan(nboot,length(InclUnits));
                resnormbootoff = nan(nboot,length(InclUnits));
                parfor bootid=1:nboot
                    bootid
                    % Evoked RFs
                    % EVoked visual responses
                    tmp = SpikeCountPerTP(InclTP,InclUnits);
                    tmp = tmp(datasample(1:size(tmp,1),size(tmp,1),'replace',true),:); %Shuffle time points!
                    RFON = arrayfun(@(X) nanmean(tmp(NewStimVecON(InclTP,X)==1,:),1),1:newY*newX,'UniformOutput',0);
                    RFON = cat(1,RFON{:});
                    RFOFF = arrayfun(@(X) nanmean(tmp(NewStimVecOFF(InclTP,X)==1,:),1),1:newY*newX,'UniformOutput',0);
                    RFOFF = cat(1,RFOFF{:});
                    
                    % Again, Z-score
                    RFON = (RFON - nanmean(RFON,1))./nanstd(RFON,[],1);
                    RFOFF = (RFOFF - nanmean(RFOFF,1))./nanstd(RFOFF,[],1);
                    
                    
                    %Smooth
                    %                     RFAll = arrayfun(@(X) smooth2a(reshape(RFAll(:,X),newY,newX),2),1:size(RFAll,2),'UniformOutput',0);
                    %                     RFAll = cat(3,RFAll{:});
                    RFON = arrayfun(@(X) smooth2a(reshape(RFON(:,X),newY,newX),2),1:size(RFON,2),'UniformOutput',0);
                    RFON = cat(3,RFON{:});
                    RFOFF = arrayfun(@(X) smooth2a(reshape(RFOFF(:,X),newY,newX),2),1:size(RFOFF,2),'UniformOutput',0);
                    RFOFF = cat(3,RFOFF{:});
                    
                    RFAll = abs(RFON) + abs(RFOFF);
                    % Interpollate nans
                    RFAll = fillmissing(RFAll,'spline');
                    RFON = fillmissing(RFON,'spline');
                    RFOFF = fillmissing(RFOFF,'spline');

                    % Z-score
                    RFAll = reshape(RFAll,newY*newX,[]);
                    RFAll = reshape((RFAll-nanmean(RFAll,1))./nanstd(RFAll,[],1),newY,newX,[]);

                    RFON = reshape(RFON,newY*newX,[]);
                    RFON = reshape((RFON-nanmean(RFON,1))./nanstd(RFON,[],1),newY,newX,[]);

                    RFOFF = reshape(RFOFF,newY*newX,[]);
                    RFOFF = reshape((RFOFF-nanmean(RFOFF,1))./nanstd(RFOFF,[],1),newY,newX,[]);
            
                    % Fit best parameters to all clusters
                    x0 = [1,0,15,0,15,0]; % Initial guess [Amplitude, x0, sigmax, y0, sigmay, angel(rad)]
                    x0 = repmat(x0,[length(InclUnits),1]);
                    [x0(:,1),id] = nanmax(reshape(RFAll,[],length(InclUnits)),[],1); % Initial guess [Amplitude, x0, sigmax, y0, sigmay, angel(rad)]
                    x0(:,2) = X(id);
                    x0(:,4) = Y(id);
                    try
                        %All
                        [thisparamsboot,restmp,~,~] = arrayfun(@(Z) lsqcurvefit(@D2GaussFunctionRot,x0(Z,:),cat(3,X,Y),RFAll(:,:,Z),lb,ub,options),...
                            1:length(InclUnits),'UniformOutput',0);
                        resnormbootall(bootid,:) = cat(1,restmp{:});

                        thisparamsboot = cat(1,thisparamsboot{:});

                        
                        % OFF
                        x0 = thisparamsboot; %replace initial parameters to save time
                        [thisparamsboot,restmp,~,~] = arrayfun(@(Z) lsqcurvefit(@D2GaussFunctionRot,x0(Z,:),cat(3,X,Y),RFON(:,:,Z),lb,ub,options),...
                            1:length(InclUnits),'UniformOutput',0);
                        resnormbooton(bootid,:)  = cat(1,restmp{:});
                        % OFF
                        [thisparamsboot,restmp,~,~] = arrayfun(@(Z) lsqcurvefit(@D2GaussFunctionRot,x0(Z,:),cat(3,X,Y),RFOFF(:,:,Z),lb,ub,options),...
                            1:length(InclUnits),'UniformOutput',0);
                        resnormbootoff(bootid,:) = cat(1,restmp{:});
                                         
                    catch ME
                        disp(ME)
                    end
                end
                
                %% Calculate p-value
                pvalsAll = arrayfun(@(X) invprctile(resnormbootall(:,X),resnormAll(X))./100,1:length(InclUnits),'UniformOutput',0); %pvalue
                pvalsAll = cat(1,pvalsAll{:});
                pvalsON = arrayfun(@(X) invprctile(resnormbooton(:,X),resnormON(X))./100,1:length(InclUnits),'UniformOutput',0); %pvalue
                pvalsON = cat(1,pvalsON{:});
                pvalsOFF = arrayfun(@(X) invprctile(resnormbootoff(:,X),resnormOFF(X))./100,1:length(InclUnits),'UniformOutput',0); %pvalue
                pvalsOFF = cat(1,pvalsOFF{:});
                %% Draw RFs for significant units only to save time
                sigunits =find(pvalsAll<0.1|pvalsON<0.1|pvalsOFF<0.1)';%[1:length(pvalsAll)]';%
                if drawthis
                    for unitid=sigunits
                        figure('name',[areashere{unitid} 'RF, Cell' num2str(cluster_id(Good_IDx(InclUnits(unitid)))) ', RFfit p = ' num2str(pvalsON(unitid))])
                        a1 = subplot(2,2,1); h=imagesc(xwidth,yheight,RFON(:,:,unitid)); set(gca,'YDir','normal'); title(['ON, p=' num2str(pvalsON(unitid))])
                        hold on;
                        ellipse(thisparamsON(unitid,3),thisparamsON(unitid,5),-thisparamsON(unitid,6),thisparamsON(unitid,2),thisparamsON(unitid,4),'r')
                        colormap gray
                        
                        a2 = subplot(2,2,2); h=imagesc(xwidth,yheight,RFOFF(:,:,unitid)); set(gca,'YDir','normal'); title(['OFF, p=' num2str(pvalsOFF(unitid))])
                        hold on;
                        ellipse(thisparamsOFF(unitid,3),thisparamsOFF(unitid,5),-thisparamsOFF(unitid,6),thisparamsOFF(unitid,2),thisparamsOFF(unitid,4),'b')
                        
                        a3 = subplot(2,2,3); h=imagesc(xwidth,yheight,RFAll(:,:,unitid)); set(gca,'YDir','normal'); title(['ON+OFF p=' num2str(pvalsAll(unitid))])
                        hold on;
                        ellipse(thisparamsAll(unitid,3),thisparamsAll(unitid,5),-thisparamsAll(unitid,6),thisparamsAll(unitid,2),thisparamsAll(unitid,4),'g')
                        
                        colormap gray
                        drawnow
                        saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,['Unit' num2str(cluster_id(Good_IDx(InclUnits(unitid)))) '_' areashere{unitid} '-RF.fig']))
                        saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,['Unit' num2str(cluster_id(Good_IDx(InclUnits(unitid)))) '_' areashere{unitid} '-RF.bmp']))
                    end
                end
                
                
                %% Organize areas
                if histoflag
                    depths = Depth2AreaPerUnit.Depth(ismember(Depth2AreaPerUnit.Cluster_ID,cluster_id(Good_IDx(InclUnits))));
                    areatmp = Depth2AreaPerUnit.Area(ismember(Depth2AreaPerUnit.Cluster_ID,cluster_id(Good_IDx(InclUnits))));
                    colperunit = cellfun(@(X) hex2rgb(X),Depth2AreaPerUnit.Color(ismember(Depth2AreaPerUnit.Cluster_ID,cluster_id(Good_IDx(InclUnits)))),'UniformOutput',0);
                    colperunit = cat(1,colperunit{:});
                else
                    areatmp = areashere;
                    depths = depth(Good_IDx(InclUnits));
                    Depth2AreaPerUnit.Depth = depths;
                    Depth2AreaPerUnit.Area = areatmp;
                end
                % Sort on depth
                [depths, sortidx] = sort(depths);
                areatmp = areatmp(sortidx);
                areatmp = strrep(areatmp,'/','');
                
                if histoflag
                    colperunit = colperunit(sortidx,:);
                else
                    colperunit = copper(length(depths));
                end
                
                
                areaopt = unique(areatmp,'stable');
                uniquearean = length(areaopt);
                % Use Allen Brain Atlas to find all areas in the recordings
                % Add units to area specific cell
                % areatmp = areatmp(~cellfun(@isempty,areatmp));
                if histoflag
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
                else
                    newareaabrev = areatmp;
                end
                
                %% Overview Figure
                col4this = copper(length(InclUnits));
                yticks = [];
                UniqueAreas = unique(newareaabrev,'stable');
                for areaid=1:length(UniqueAreas)
                    yticks = [yticks find(ismember(newareaabrev,UniqueAreas{areaid})==1,1,'first')];
                end
                [yticks,sortid]=sort(yticks,'ascend');
                UniqueAreas = UniqueAreas(sortid);
                %                        pvalsON = nan(1,length(InclUnits));
                %                     pvalsOFF = nan(1,length(InclUnits));
                %Normalize pvalcorr
                figure('name','Overview RFs')
                subplot(4,2,[1:2])
                for unitid=1:length(InclUnits)
                    
                    if pvalsON(unitid)>0.1
                        continue
                    end
                    h=ellipse(thisparamsON(unitid,3),thisparamsON(unitid,5),-thisparamsON(unitid,6),thisparamsON(unitid,2),thisparamsON(unitid,4));
                    h.Color = col4this(unitid,:);
                    h.LineWidth=2;
                    hold on
                end
                xlim([-200 200])
                ylim([-50 50])
                box off
                title('RF ON distribution')
                
                subplot(4,2,[3:4])
                for unitid=1:length(InclUnits)
                    
                    if pvalsOFF(unitid)>0.1
                        continue
                    end
                    h=ellipse(thisparamsOFF(unitid,3),thisparamsOFF(unitid,5),-thisparamsOFF(unitid,6),thisparamsOFF(unitid,2),thisparamsOFF(unitid,4));
                    h.Color = col4this(unitid,:);
                    h.LineWidth=2;
                    hold on
                end
                xlim([-200 200])
                ylim([-50 50])
                box off
                title('RF OFF distribution')
                
                subplot(4,2,[5:6])
                for unitid=1:length(InclUnits)
                    
                    if pvalsAll(unitid)>0.1
                        continue
                    end
                    h=ellipse(thisparamsAll(unitid,3),thisparamsAll(unitid,5),-thisparamsAll(unitid,6),thisparamsAll(unitid,2),thisparamsAll(unitid,4));
                    h.Color = col4this(unitid,:);
                    h.LineWidth=2;
                    hold on
                end
                xlim([-200 200])
                ylim([-50 50])
                box off
                title('RF All distribution')
                
                subplot(3,2,[5,6])
                imagesc([1:length(InclUnits)],1,[1:length(InclUnits)])
                colormap(col4this)
                %                     imagesc(uniquedepths,1,uniquedepths')
                set(gca,'XTick',(yticks),'XTickLabels',UniqueAreas,'XTickLabelRotation',25,'YTickLabel','')
                
                %                     set(gca,'XTick',uniquedepths(yticks),'XTickLabels',UniqueAreas,'XTickLabelRotation',25,'YTickLabel','')
                box off
                makepretty
                %                     colormap(jet)
                saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,['RFOverview.fig']))
                saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,['RFOverview.bmp']))
                
                
                %% Save data
                RFInfo = table;
                RFInfo.ClusID = cluster_id(Good_IDx(InclUnits));
                RFInfo.paramsON = thisparamsON;
                RFInfo.paramsOFF = thisparamsOFF;
                RFInfo.paramsAll = thisparamsAll;
                
                RFInfo.PvalON = pvalsON;
                RFInfo.PvalOFF = pvalsOFF;
                RFInfo.PvalsAll = pvalsAll;
                RFInfo.RFON = permute(RFON,[3,1,2]);
                RFInfo.RFOFF = permute(RFOFF,[3,1,2]);
                RFInfo.RFAll = permute(RFAll,[3,1,2]);
                RFInfo.Depth = depth(Good_IDx(InclUnits));
                RFInfo.Shank = Shank(Good_IDx(InclUnits));
                
                
                ExpInfo.azimuth = xwidth;
                ExpInfo.elevation = yheight;
                ExpInfo.Protocol = Protocol;
                
                
                save(fullfile(SaveRFDir,MiceOpt{midx},thisdate,thisses,thisprobe,'RFData.mat'),'RFInfo','ExpInfo','Depth2AreaPerUnit','clusinfo','SpikeCountPerTP','-v7.3')
                
                clear RFInfo
                clear ExpInfo
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
                clear Depth2AreaPerUnit
            end
        end
    end
    close all
    %% Load all results to make big overview per mouse
    AllRFDat=[];
    for didx = 1:length(Dates4Mouse)
        % Within folders, look for 'RF mapping sessions'
        AllRFDatFiles = dir(fullfile(SaveRFDir,MiceOpt{midx},'*','*','*','RFData.mat'));
        for idx=1:length(AllRFDatFiles)
            % Load Data
            tmp = load(fullfile(AllRFDatFiles(idx).folder,AllRFDatFiles(idx).name));
            
            % Save in big struct
            if isempty(AllRFDat)
                
                AllRFDat = tmp.RFInfo;
                AllRFDat = tall(AllRFDat); %Make tall to prevent memory issues
                
            else
                % Concatenate
                AllRFDat = cat(1,AllRFDat,tmp.RFInfo);
                
            end
            
        end
    end
    
    %% Plot
    depth = gather(AllRFDat.Depth);
    Shank = gather(AllRFDat.Shank);
    paramsAll = gather(AllRFDat.paramsAll);
    pvalAll = gather(AllRFDat.PvalsAll);
    [depth,unitsort] = sort(depth,'ascend');
    col4this = copper(length(depth));
    yticks = [];
    newareaabrev = arrayfun(@(X) [num2str(round(depth(X)./100)*100) 'um'],1:length(depth),'UniformOutput',0);
    
    UniqueAreas = unique(newareaabrev,'stable');
    for areaid=1:length(UniqueAreas)
        yticks = [yticks find(ismember(newareaabrev,UniqueAreas{areaid})==1,1,'first')];
    end
    [yticks,sortid]=sort(yticks,'ascend');
    UniqueAreas = UniqueAreas(sortid);
    %                        pvalsON = nan(1,length(InclUnits));
    %                     pvalsOFF = nan(1,length(InclUnits));
    %Normalize pvalcorr
    figure('name','Overview RFs')
    subplot(4,2,[1:6])
    for unitid=1:length(depth)
        if pvalAll(unitsort(unitid))>0.05
            continue
        end
        
        h=ellipse(paramsAll(unitsort(unitid),3),paramsAll(unitsort(unitid),5),-paramsAll(unitsort(unitid),6),paramsAll(unitsort(unitid),2),paramsAll(unitsort(unitid),4));
        h.Color = col4this(unitid,:);
        h.LineWidth=2;
        hold on
    end
    xlim([-200 200])
    ylim([-50 50])
    box off
    title('RF All distribution')
    
    
    subplot(3,2,[5,6])
    imagesc([1:length(depth)],1,[1:length(depth)])
    colormap(col4this)
    %                     imagesc(uniquedepths,1,uniquedepths')
    set(gca,'XTick',(yticks),'XTickLabels',UniqueAreas,'XTickLabelRotation',25,'YTickLabel','')
    
    %                     set(gca,'XTick',uniquedepths(yticks),'XTickLabels',UniqueAreas,'XTickLabelRotation',25,'YTickLabel','')
    box off
    makepretty
    %                     colormap(jet)
    saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},['RFOverview.fig']))
    saveas(gcf,fullfile(SaveRFDir,MiceOpt{midx},['RFOverview.bmp']))
    
    
    
    
end

