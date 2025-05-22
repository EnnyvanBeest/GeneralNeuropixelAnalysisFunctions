if ~exist('Actualtime')
    loadTimeline
end

%Find directory raw data
if exist('Params') && isfield(Params,'RawDataPaths') && ~isempty(Params.RawDataPaths)
    if length(Params.RawDataPaths) > 1
        if isstruct(Params.RawDataPaths{1})
            Params.RawDataPaths = cellfun(@(X) fullfile(X.folder,X.name),Params.RawDataPaths,'Uni',0);
        end
        Params.RawDataPaths = unique(Params.RawDataPaths,'stable');
        lfpDAll = cell2mat(cellfun(@(X) dir(X),Params.RawDataPaths,'Uni',0));

    elseif isstruct(Params.RawDataPaths{1})
        lfpDAll = Params.RawDataPaths{1};
    else
        lfpDAll = dir(Params.RawDataPaths{1});
    end
else
    try %Reading from params file in kilosort directory
        lfpDAll = cell(1,numel(myKsDir));
        for ksid = 1:numel(myKsDir)
            spikeStruct = loadParamsPy(fullfile(myKsDir{ksid}, 'params.py'));
            rawD = spikeStruct.dat_path;
            rawD = strsplit(rawD, ',');
            for rid = 1:length(rawD)
                rawD{rid} = rawD{rid}(strfind(rawD{rid}, '"') + 1:end);
                rawD{rid} = rawD{rid}(1:strfind(rawD{rid}, '"') - 1);
                rawD{rid} = dir(rawD{rid});
                if isempty(rawD{rid})
                    rawD{rid} = dir(strrep(rawD{rid}, 'bin', 'cbin'));
                end
            end
            rawD = cat(2, rawD{:});
            lfpDAll{ksid} = rawD;
        end
        lfpDAll = [lfpDAll{:}];
    catch


        lfpDAll = dir(fullfile(myLFDir,'*','*',['*imec' num2str(probeid-1) '.ap.*bin'])); % LFP file from spikeGLX specifically
    end
end

clear syncDatImec
syncchanmissing = 1;
for recordingsessionidx = 1:length(lfpDAll)
    lfpD = lfpDAll(recordingsessionidx);

    spikeTimestmp = st(RecSes==recordingsessionidx);            % Spike times according to IMEC

    spikeTimesCorrected = nan(size(spikeTimestmp)); % corrected spike times (to timeline)
    spikeTimesAlignedRecording = nan(size(st));

    warningflag = 0;

    abortthissession = 0;
    AllDelayIdx = [];


    %Sort sessions
    allsess = dir(fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},thisdate));
    sesnrs = (cellfun(@(X) str2num(X),{allsess(:).name},'UniformOutput',0));
    sesnrs(cellfun(@isempty,sesnrs))={[1000]};
    [~,idx] = sort([sesnrs{:}]);
    allsess = allsess(idx);
    % NIDQ FILE
    if nidq_sync_used(midx)
        tmpfile = strsplit(lfpDAll.folder,filesep);
        tmpfile = fullfile('\\',tmpfile{2:end-2});

        mySyncFile = dir(fullfile(tmpfile,'**','*nidq.bin'))
        flag = 0;
        idx=1;
        while ~flag
            if idx>length(mySyncFile)
                warning('No matching data, skip')
                break
            end
            if ~isempty(strfind(mySyncFile(idx).name,'Concat'))
                idx = idx+1;
                continue
            end

            [dataArray,samplerate,stime,endtime] = ReadSGLXData(mySyncFile(idx).name,mySyncFile(idx).folder,length(Actualtime)./tmSR,find(ismember({allsess(:).name},thisses))); %1 To sync with NP, 2 acquire live, 3 flipper

            % Match this to the flipper signal
            FlipperGLX = dataArray(3,:);
            FlipperGLX = downsample(FlipperGLX,round(samplerate/tmSR)); %Make digital
            FlipperGLX(FlipperGLX<2)=0;
            FlipperGLX(FlipperGLX>2)=1;

            indx=1:length(Flippertimeline);
            % Initial estimate of delay
            delayestimate = finddelay(Flippertimeline(indx),FlipperGLX);
            indx2 = indx+delayestimate;
            indx(indx2<1)=[];
            indx2(indx2<1)=[];
            indx(indx2>length(FlipperGLX))=[];
            indx2(indx2>length(FlipperGLX))=[];
            if corr(Flippertimeline(indx),FlipperGLX(indx2)')<0.85 %Test if signals correlate enough
                idx=idx+1;
                continue
            end

            syncchanmissing = 0;

            % NIDQ PULSE (=same clock as FlipperGLX)
            %                 mySyncFile = dir(fullfile(myLFDir,'*','*nidq.bin'));
            syncChanIndex = 1; %NPSync
            syncDatNIDQ = dataArray(syncChanIndex,:);%extractSyncChannel(mySyncFile(1).folder, 3, syncChanIndex); % Further digitize
            tmpmean = nanmean(syncDatNIDQ(:));
            syncDatNIDQ(syncDatNIDQ<tmpmean)=0;
            syncDatNIDQ(syncDatNIDQ>tmpmean)=1;
            %                 [nidqmeta] = (ReadMeta2(mySyncFile(1).folder));
            %Align to timeline
            syncDatNIDQ = downsample(syncDatNIDQ,round(samplerate/tmSR));
            %                                 syncDatNIDQTime = [1:length(syncDatNIDQ)]./str2double(nidqmeta.niSampRate);

            %IMEC PULSE
            syncChanIndex = 385;
            nChansInFile = 384;
            syncDatImec = extractSyncChannel(lfpD.folder, nChansInFile, syncChanIndex); %Make sure extract sync channel read in ap data, not lf
            [Imecmeta] = (ReadMeta2(lfpD.folder));

            % normalize between 0 and 1 (it's now in uint16 till 64)
            tmpmean = nanmean(syncDatImec(:));
            syncDatImec(syncDatImec<tmpmean)=0;
            syncDatImec(syncDatImec>tmpmean)=1;

            startidx = floor(stime*str2double(Imecmeta.imSampRate));
            if startidx<1
                startidx=1;
            end
            endidx = ceil(endtime*str2double(Imecmeta.imSampRate));
            if endidx>length(syncDatImec)
                endidx=length(syncDatImec)
            end
            syncDatImec = syncDatImec(startidx:endidx);
            %Align to timeline
            syncDatImec = downsample(syncDatImec,round(str2double(Imecmeta.imSampRate)/tmSR));

            indx(indx2>length(syncDatImec))=[];
            indx2(indx2>length(syncDatImec))=[];
            % Initial estimate of delay
            delayestimate2 = finddelay(syncDatNIDQ(indx2),syncDatImec(indx2));
            indx3 = indx2+delayestimate2;
            indx(indx3<1) = [];
            indx2(indx3<1)=[];
            indx3(indx3<1)=[];
            indx(indx3>length(syncDatImec))=[];
            indx2(indx3>length(syncDatImec))=[];
            indx3(indx3>length(syncDatImec))=[];
            if corr(double(syncDatNIDQ(indx2)'),double(syncDatImec(indx3)'))<0.9
                idx=idx+1;
                continue
            end
            flag=1;
        end
        if flag
            figure;
            per2check = 5;
            % align timelines every x seconds (align GLX to Timeline)
            for i = 1:per2check:max(ceil(Actualtime))
                indx=find(Actualtime>=i-1 & Actualtime<i+per2check-1);
                indx(indx>length(Actualtime))=[];
                indx(indx>length(FlipperGLX))=[];

                if isempty(indx)
                    break
                end

                if exist('h1')
                    delete(h1)
                    delete(h2)

                end

                if exist('h3')
                    delete(h3)
                end
                if exist('h4')
                    delete(h4)
                end
                % Check every x timepoints what the difference is between the
                % two flippers
                indx2 = indx+delayestimate;
                indx(indx2<1)=[];
                indx2(indx2<1)=[];
                indx(indx2>length(FlipperGLX))=[];
                indx2(indx2>length(FlipperGLX))=[];
                if isempty(indx2)
                    continue
                end
                delayindx = finddelay(Flippertimeline(indx),FlipperGLX(indx2));

                %D = FINDDELAY(X,Y), where X and Y are row or column vectors of length
                %   LX and LY, respectively, returns an estimate of the delay D between X
                %   and Y, where X serves as the reference vector. If Y is delayed with
                %   respect to X, D is positive. If Y is advanced with respect to X, D is
                %   negative.

                indx2 = indx2+delayindx;
                indx(indx2<1)=[];
                indx2(indx2<1)=[];
                indx(indx2>length(FlipperGLX))=[];
                indx2(indx2>length(FlipperGLX))=[];
                if isempty(indx2)
                    continue
                end
                h1 = plot(Actualtime(indx),Flippertimeline(indx),'b-');
                hold on
                h2 = plot(Actualtime(indx),FlipperGLX(indx2)+1,'r-');


                indx3 = indx2+delayestimate2;
                indx(indx3<1)=[];
                indx2(indx3<1)=[];
                indx3(indx3<1)=[];
                indx(indx3>length(syncDatImec))=[];
                indx2(indx3>length(syncDatImec))=[];
                indx3(indx3>length(syncDatImec))=[];
                % sync channel data (ch. 385 NP) to Nidq channel 1
                delayindx2 =finddelay(syncDatNIDQ(indx2),syncDatImec(indx3));

                % +delayindx2
                indx3 = indx3+delayindx2;
                indx(indx3<1)=[];
                indx2(indx3<1)=[];
                indx3(indx3<1)=[];
                indx(indx3>length(syncDatImec))=[];
                indx2(indx3>length(syncDatImec))=[];
                indx3(indx3>length(syncDatImec))=[];

                h3 = plot(Actualtime(indx),syncDatNIDQ(indx2)+2,'r-');
                h4 = plot(Actualtime(indx),syncDatImec(indx3)+3,'k-');

                %sanity check
                %                     figure;
                %                     indx=abs((-delayindx2-delayestimate2-delayindx-delayestimate))+1:abs((-delayindx2-delayestimate2-delayindx-delayestimate))+5000;
                %                     plot(syncDatImec(indx),'k-')
                %                     hold on
                %                     plot(syncDatNIDQ(indx-delayindx2-delayestimate2)+1,'r-')
                %                     plot(FlipperGLX(indx-delayindx2-delayestimate2)+2,'r-')
                %                     plot(Flippertimeline(indx-delayindx2-delayestimate2-delayindx-delayestimate)+3,'b-')
                %
                %Convert back to time in imec data
                TL2ImecTime = (delayindx2+delayestimate2+delayindx+delayestimate)./tmSR;
                Imec2TLTime = (-delayindx2-delayestimate2-delayindx-delayestimate)./tmSR;

                %Find spikes in IMEC time, that fall in this window in
                %TIMELINE space
                spikeindx = find(spikeTimestmp>=TL2ImecTime+stime+(i-1)&spikeTimestmp<=TL2ImecTime+stime+(i+per2check-1));
                if ~isempty(spikeindx)
                    spikeTimesCorrected(spikeindx)=spikeTimestmp(spikeindx)-stime+Imec2TLTime;
                end
                drawnow
                AllDelayIdx = [AllDelayIdx -stime+Imec2TLTime];
            end

        else
            disp('Cannot find correlating data')
            abortthissession=1;
        end
        % Just once align all
        if ~syncchanmissing

            spikeTimesAlignedRecording = spikeTimestmp+nanmean(AllDelayIdx);
            if exist('syncDatImec')
                syncDatImec = circshift(syncDatImec,round(nanmean(AllDelayIdx).*tmSR));
            end
        end

    else % Flipper directly in sync channel;

        indx=1:length(Flippertimeline);

        %IMEC PULSE
        [Imecmeta] = ReadMeta2(lfpD.folder,'ap');
        nchan = strsplit(Imecmeta.acqApLfSy,',');
        nChansInFile = str2num(nchan{1})+str2num(nchan{3});
        syncDatImec = extractSyncChannel(lfpD.folder, nChansInFile, nChansInFile); %Last channel is sync
        if isstruct(syncDatImec)
            syncDatImec = syncDatImec.syncdata;
        end

        % normalize between 0 and 1 (it's now in uint16 till 64)
        tmpmean = nanmean(syncDatImec(:));
        if tmpmean>0
            syncDatImec(syncDatImec<tmpmean)=0;
            syncDatImec(syncDatImec>tmpmean)=1;
        else
            syncDatImec(syncDatImec>tmpmean)=1;
            syncDatImec(syncDatImec<tmpmean)=0;
        end

        %Align to timeline
        if isempty(syncDatImec)
            keyboard
            continue
        end
        firstval = syncDatImec(1);

        syncDatImec = downsample(syncDatImec,round(str2double(Imecmeta.imSampRate)/tmSR));

        % Aproximate start of session relative to ephys
        %             fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,'SpikeData.mat')
        allsess = dir(fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},thisdate));
        sesnrs = (cellfun(@(X) str2num(X),{allsess(:).name},'UniformOutput',0));
        sesnrs(cellfun(@isempty,sesnrs))={[1000]};
        [~,idx] = sort([sesnrs{:}]);
        allsess = allsess(idx);
        totaldur = 0;

        % now find in syncDatImec when new session started (first 1 after period of 0's)
        WindowSize = round(5.*tmSR); %Should be nothing for at least 5 seconds
        filtereddat = filter((1/tmSR)*ones(1,WindowSize),1,syncDatImec);
        filtereddat(filtereddat==max(filtereddat))=0;

        okay = 0;
        while ~okay
            probablestartpoint = find(filtereddat(round(totaldur*tmSR)+1:end)<0.1,1,'first')+totaldur*tmSR;
            probablestartpoint = find(filtereddat(probablestartpoint:end)>0.1,1,'first')+probablestartpoint;

            if isempty(probablestartpoint) || probablestartpoint>length(filtereddat)
                warningflag = 1;
                okay = 1;
                delayestimate = 0;
                continue
            end

            %         if ~exist('probablestartpoint','var') || isempty(probablestartpoint)
            %             probablestartpoint = totaldur*tmSR;
            %         end
            indx2 = indx+round(probablestartpoint);
            indx(indx2<1) = [];
            indx2(indx2<1)=[];
            indx(indx2>length(syncDatImec))=[];
            indx2(indx2>length(syncDatImec))=[];
            % Initial estimate of delay
            delayestimate = finddelay(Flippertimeline(single(indx)),syncDatImec(single(indx2)))+probablestartpoint;
            indx2 = indx+delayestimate;
            syncchanmissing = 0;

            indx(indx2<1) = [];
            indx2(indx2<1)=[];
            indx(indx2>length(syncDatImec))=[];
            indx2(indx2>length(syncDatImec))=[];
            if size(syncDatImec,1)<size(syncDatImec,2)
                syncDatImec = syncDatImec';
            end
            if corr(double(Flippertimeline(single(indx))),double(syncDatImec(single(indx2))))<0.75 || isnan(corr(double(Flippertimeline(single(indx))),double(syncDatImec(single(indx2)))))
                tmpfig = figure; plot(Flippertimeline(single(indx))); hold on; plot(syncDatImec(single(indx2))+1)
                if probablestartpoint<length(filtereddat)
                    totaldur = probablestartpoint/tmSR+5;%add a few seconds
                    indx=1:length(Flippertimeline);
                else
                    warningflag = 1;
                    okay = 1;
                end
            else
                okay=1;
            end
        end

        figure;

        per2check = 5;
        % align timelines every x seconds (align GLX to Timeline)
        for i = 1:per2check:max(ceil(Actualtime))
            if syncchanmissing
                break
            end
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
            delayindx = finddelay(Flippertimeline(single(indx)),syncDatImec(single(indx2)));

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
            h1 = plot(Actualtime(single(indx)),Flippertimeline(single(indx)),'b-');
            hold on
            h2 = plot(Actualtime(single(indx)),syncDatImec(single(indx2))+1,'r-');

            if corr(double(Flippertimeline(single(indx))),double(syncDatImec(single(indx2))))>0.9
                warningflag=0;
            elseif warningflag && corr(double(Flippertimeline(single(indx))),double(syncDatImec(single(indx2))))<0.7
                warningflag = 1;
            end
            %sanity check
            if warningflag && i > 5
                figure;
                indx=abs((-delayindx-delayestimate))+1:abs((-delayindx-delayestimate))+5000;
                plot(syncDatImec(single(indx)),'k-')
                hold on
                plot(Flippertimeline(single(indx-delayindx-delayestimate))+1,'b-')
                drawnow
                AcceptableYN = '';
                while ~ismember(AcceptableYN,{'y','n'})
                    AcceptableYN = 'n';%input('Is this alignment acceptable? (y/n)','s');
                    if strcmpi(AcceptableYN,'n')
                        disp('Not acceptable, see if alignment can be done with spike data')
                        syncchanmissing = 1;
                        break
                    elseif strcmp(AcceptableYN,'y')
                        disp('Acceptable, continuing...')
                    end
                end
            end
            if syncchanmissing
                break
            end

            %Convert back to time in imec data
            TL2ImecTime = (delayindx+delayestimate)./tmSR;
            Imec2TLTime = -TL2ImecTime;
            %Sanity check
            %         figure;
            %         indx=abs((-delayindx-delayestimate))+1:abs((-delayindx-delayestimate))+5000;
            %         plot(syncDatImec(indx),'k-')
            %         hold on
            %         plot(Flippertimeline(indx-delayindx-delayestimate)+1,'b-')
            %
            %         Find spikes in IMEC time, that fall in this window in
            %TIMELINE space
            spikeindx = find(spikeTimestmp>=TL2ImecTime+(i-1)&spikeTimestmp<=TL2ImecTime+(i+per2check-1));
            if ~isempty(spikeindx)
                if any(spikeindx>numel(spikeTimestmp))
                    keyboard
                end
                spikeTimesCorrected(spikeindx)=spikeTimestmp(spikeindx)+Imec2TLTime;
            end
            AllDelayIdx = [AllDelayIdx Imec2TLTime];

            drawnow
        end
        % Just once align all
        if ~syncchanmissing
            spikeTimesAlignedRecording = spikeTimestmp+nanmean(AllDelayIdx);
            if exist('syncDatImec') && ~isempty(AllDelayIdx)
                syncDatImec = circshift(syncDatImec,round(nanmean(AllDelayIdx).*tmSR));
            end
        end
    end


    if ~syncchanmissing
        break
    end
end

if ~syncchanmissing

tmp = detrend(double(syncDatImec));
tmp = (tmp-nanmin(tmp(:)))./(nanmax(tmp(:))-nanmin(tmp(:)));
tresh=0.5;
tmp(tmp>=tresh)=1;
tmp(tmp<tresh) = 0;

% Convert to 1 for each flip (so dark and light square change is 1)
tmpdiff = [0; abs(diff(tmp(:)))];
tresh = 0.05;
tmpdiff(tmpdiff>=tresh)=1;
tmpdiff(tmpdiff~=1)=0;

% If ITI times are unknown, take at 20 seconds
ITITimes = 20;

% sampling rate of time line
WindowSize = round(0.8*(nanmin(ITITimes).*tmSR));

% Upperlimit for refresh rate frames measured
refreshrateInFrames = quantile((diff(find(tmpdiff==1))),0.5);

% Sliding window with average ITI times to find them approximately
filtereddat = filter((1/WindowSize)*ones(1,WindowSize),1,tmpdiff);
diffTime = [nan; diff(filtereddat(:))];

% figure; ax1 = subplot(3,1,1); plot(tmpdiff); title('syncDatImec'); hold on;
% ax2 = subplot(3,1,2); plot(filtereddat); title('Filtered syncDatImec'); hold on
% ax3 = subplot(3,1,3); plot(diffTime); title('Change in Filtered syncDatImec'); hold on
% linkaxes([ax1,ax2,ax3],'x')

tmpdif = [0; diff(filtereddat(:)==0)];
SessionStart = find(tmpdif==-1);
SessionStop = cell2mat(arrayfun(@(X) find(filtereddat(X:end)==0,1,'first')+X-1-WindowSize,SessionStart,'UniformOutput',0));

end