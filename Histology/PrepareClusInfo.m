% Create saving directory
if ~exist(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe))
    mkdir(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe))
end

% Check for multiple subfolders?
subsesopt = dir(myKsDir);
subsesopt(~[subsesopt.isdir])=[];

%% Initialize everything
sp = cell(1,0);
channelmap=[];
channelpos = [];
cluster_id = [];
Label = [];
Good_ID = [];
depth = [];
channel = [];
Shank=[];
recses = [];
countid=1;
% figure;
cols = jet(length(subsesopt));
for subsesid=1:length(subsesopt)
    if isempty(dir(fullfile(subsesopt(subsesid).folder,subsesopt(subsesid).name,'*.npy')))
        continue
    end
    
    thissubses = str2num(subsesopt(subsesid).name)
    %% Load Spike Data
    sp{countid} = loadKSdir(fullfile(subsesopt(subsesid).folder,subsesopt(subsesid).name),params); % Load Spikes with PCs
    [sp{countid}.spikeAmps, sp{countid}.spikeDepths, sp{countid}.templateDepths, sp{countid}.templateXpos, sp{countid}.tempAmps, sp{countid}.tempsUnW, sp{countid}.templateDuration, sp{countid}.waveforms] = templatePositionsAmplitudes(sp{countid}.temps, sp{countid}.winv, sp{countid}.ycoords, sp{countid}.xcoords, sp{countid}.spikeTemplates, sp{countid}.tempScalingAmps); %from the spikes toolbox
    
    
    %% Channel data
    myClusFile = dir(fullfile(subsesopt(subsesid).folder,subsesopt(subsesid).name,'channel_map.npy'));
    channelmaptmp = readNPY(fullfile(myClusFile(1).folder,myClusFile(1).name));
    
    myClusFile = dir(fullfile(subsesopt(subsesid).folder,subsesopt(subsesid).name,'channel_positions.npy'));
    channelpostmp = readNPY(fullfile(myClusFile(1).folder,myClusFile(1).name));
    if length(channelmaptmp)<length(channelpostmp)
        channelmaptmp(end+1:length(channelpostmp))=length(channelmaptmp):length(channelpostmp)-1;
    end
    %% Is it correct channelpos though...?    
    myLFDir = fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},thisdate,'ephys');
    lfpD = dir(fullfile([myLFDir '*'], '**\*.lf.*bin')); % ap file from spikeGLX specifically
    if isempty(lfpD)
    lfpD = dir(fullfile([myLFDir '*'], '**\*.ap.*bin')); % ap file from spikeGLX specifically
    end
    
    if isempty(lfpD)
        disp('No LFP data found')
    elseif length(lfpD)>length(subksdirs)
        disp('Just take data from the last recording')
        lfpD = lfpD((thissubses));
    elseif length(lfpD)<length(subksdirs)
        disp('Should be a different amount of probes?')
        keyboard
    else
        lfpD = lfpD(probeid);
    end
    
    if strcmp(ProbeType,'2_4S')
        channelpostmp = ChannelIMROConversion(lfpD.folder,1); % For conversion when not automatically done
    end
    %% Load Cluster Info
    myClusFile = dir(fullfile(subsesopt(subsesid).folder,subsesopt(subsesid).name,'cluster_info.tsv'));
    if isempty(myClusFile)
        disp('This data is not yet curated with phy!!')
        curratedflag=0;
        myClusFile = dir(fullfile(subsesopt(subsesid).folder,subsesopt(subsesid).name,'cluster_group.tsv'));
        clusinfo = tdfread(fullfile(myClusFile(1).folder,myClusFile(1).name));
        cluster_id = cat(1,cluster_id,clusinfo.cluster_id);
        tmpLabel = char(length(clusinfo.cluster_id));
        KSLabelfile = tdfread(fullfile(subsesopt(subsesid).folder,subsesopt(subsesid).name,'cluster_KSLabel.tsv'));
        tmpLabel(ismember(clusinfo.cluster_id,KSLabelfile.cluster_id)) = KSLabelfile.KSLabel(ismember(KSLabelfile.cluster_id,clusinfo.cluster_id));
        Label = [Label,tmpLabel];
        Good_ID = [Good_ID,ismember(tmpLabel,'g')]; %Identify good clusters
        
        % Find depth and channel
        depthtmp = nan(length(clusinfo.cluster_id),1);
        xtmp = nan(length(clusinfo.cluster_id),1);
        channeltmp = nan(length(clusinfo.cluster_id),1);
        for clusid=1:length(depthtmp)
            depthtmp(clusid)=round(sp{countid}.templateDepths(clusid));%round(nanmean(sp{countid}.spikeDepths(find(sp{countid}.clu==clusid-1))));
            xtmp(clusid)=sp{countid}.templateXpos(clusid);
            [~,minidx] = min(cell2mat(arrayfun(@(X) pdist(cat(1,channelpostmp(X,:),[xtmp(clusid),depthtmp(clusid)]),'euclidean'),1:size(channelpostmp,1),'UniformOutput',0)));
            try
                channeltmp(clusid) = channelmaptmp(minidx);
            catch
                channeltmp(clusid)=channeltmp(minidx-1);
            end

        end
        depth = cat(1,depth, depthtmp);
        channel = cat(1,channel,channeltmp);
        
    else
        CurationDone = 1;
        save(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,'CuratedResults.mat'),'CurationDone')
        clusinfo = tdfread(fullfile(myClusFile(1).folder,myClusFile(1).name));
        
        curratedflag=1;
        if isfield(clusinfo,'id')
            cluster_id = [cluster_id,clusinfo.id];
        elseif isfield(clusinfo,'cluster_id')
            cluster_id=[cluster_id,clusinfo.cluster_id];
        else
            keyboard
            disp('Someone thought it was nice to change the name again...')
        end
        
        % make sure cluster_id's match
        myOtherClusFile = dir(fullfile(myKsDir,'cluster_group.tsv'));%
        cluster_groupOri = tdfread(fullfile(myOtherClusFile(1).folder,myOtherClusFile(1).name));
        
        KSLabel = clusinfo.KSLabel;
        Label = [Label,clusinfo.group]; % You want the group, not the KSLABEL!
        depth = [depth,clusinfo.depth];
        channel = [channel,clusinfo.ch];
        Good_ID = [Good_ID,ismember(cellstr(clusinfo.group),'good')]; %Identify good clusters 
        channeltmp = clusinfo.ch;
    end
    
    ypostmp = channelpostmp(:,2);
    xpostmp = channelpostmp(:,1);
    xposopt = (floor(xpostmp./250));% Assuming no new shank if not at least 100 micron apart
    Shanktmp = floor(xpostmp(channeltmp+1)./250);
    
%     [~,minid] = arrayfun(@(X) (abs(floor(xpostmp(X)./250)-xposopt)),channeltmp+1,'UniformOutput',0);
    Shank = cat(1,Shank,Shanktmp); 

    recses = cat(1, recses, repmat(countid,length(channeltmp),1));
    
    channelpos = cat(1,channelpos,channelpostmp);
    channelmap = cat(1,channelmap,channelmaptmp);
    
%     scatter(xtmp,depthtmp,10,cols(countid,:))
%     hold on
%     drawnow
    countid=countid+1;

end
sp = [sp{:}];
% Find correct dataset index
spikeTimes =cat(1,sp(:).st);
spikeSites = cat(1,sp(:).spikeTemplates);
spikeCluster = cat(1,sp(:).clu);
spikeAmps = cat(1,sp(:).spikeAmps);
spikeDepths = cat(1,sp(:).spikeDepths);
spikeShank = nan(length(spikeCluster),1);
ShankOpt = unique(Shank);
for shid = 1:length(ShankOpt)
    spikeShank(ismember(spikeCluster,cluster_id(Shank==ShankOpt(shid)))) = shid;
end
templateDepths = cat(1,sp(:).templateDepths);
tempAmps = cat(1,sp(:).tempAmps);
tempsUnW = cat(1,sp(:).tempsUnW);
templateDuration = cat(1,sp(:).templateDuration);
waveforms = cat(1,sp(:).waveforms);


%  [sp{countid}.spikeAmps, sp{countid}.spikeDepths, sp{countid}.templateDepths, sp{countid}.tempAmps, sp{countid}.tempsUnW, sp{countid}.templateDuration, sp{countid}.waveforms] 
ShankOpt = unique(Shank);
ShankID = nan(size(Shank));
for id = 1:length(ShankOpt)
    ShankID(Shank==ShankOpt(id))=id;
end
clusinfo.Shank = Shank;
clusinfo.ShankID = ShankID;
Good_IDx = find(Good_ID);
clusinfo.ch = channel;
clusinfo.depth = depth;
clusinfo.cluster_id = cluster_id;
clusinfo.group = Label;
clusinfo.Good_ID = Good_ID;
