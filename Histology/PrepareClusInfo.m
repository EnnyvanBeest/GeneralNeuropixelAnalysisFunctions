% Create saving directory
if ~exist(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe))
    mkdir(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe))
end

myClusFile = dir(fullfile(myKsDir,'cluster_info.tsv'));
if isempty(myClusFile)
    disp('This data is not yet curated with phy!!')
    curratedflag=0;
    myClusFile = dir(fullfile(myKsDir,'cluster_group.tsv'));
    clusinfo = tdfread(fullfile(myClusFile(1).folder,myClusFile(1).name));
    cluster_id = clusinfo.cluster_id;
    KSLabel = char(length(cluster_id));
    KSLabelfile = tdfread(fullfile(myKsDir,'cluster_KSLabel.tsv'));
    KSLabel(ismember(cluster_id,KSLabelfile.cluster_id)) = KSLabelfile.KSLabel(ismember(KSLabelfile.cluster_id,cluster_id));
    Good_ID = ismember(cellstr(KSLabel),'good'); %Identify good clusters
    depth = nan(length(cluster_id),1);
    for clusid=1:length(depth)
        depth(clusid)=round(nanmean(spikeDepths(find(spikeCluster==clusid-1))));
    end
    myClusFile = dir(fullfile(myKsDir,'channel_map.npy'));
    channelmap = readNPY(fullfile(myClusFile(1).folder,myClusFile(1).name));
    myClusFile = dir(fullfile(myKsDir,'channel_positions.npy'));
    channelpos = readNPY(fullfile(myClusFile(1).folder,myClusFile(1).name));
    channelpos = channelpos(:,2);
    
    channel = nan(length(cluster_id),1);
    for clusid=1:length(channel)
        [minval,idx]=min(abs(depth(clusid)-channelpos));
        channel(clusid) = channelmap(idx);
    end
    clusinfo.ch = channel;
    clusinfo.depth = depth;
    clusinfo.id = cluster_id;
else
    CurationDone = 1;
    save(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,'CuratedResults.mat'),'CurationDone')
    clusinfo = tdfread(fullfile(myClusFile(1).folder,myClusFile(1).name));
    curratedflag=1;
    if isfield(clusinfo,'id')
        cluster_id = clusinfo.id;
    elseif isfield(clusinfo,'cluster_id')
        cluster_id=clusinfo.cluster_id;
    else
       keyboard
       disp('Someone thought it was nice to change the name again...')       
    end
    KSLabel = clusinfo.KSLabel;
    depth = clusinfo.depth;
    channel = clusinfo.ch;
    Good_ID = ismember(cellstr(KSLabel),'good'); %Identify good clusters
end

% Find shank
myClusFile = dir(fullfile(myKsDir,'channel_map.npy'));
channelmap = readNPY(fullfile(myClusFile(1).folder,myClusFile(1).name));

myClusFile = dir(fullfile(myKsDir,'channel_positions.npy'));
channelpos = readNPY(fullfile(myClusFile(1).folder,myClusFile(1).name));

% Shank options
xpos = channelpos(:,1);
Shank = ceil(xpos(channel+1)/100); % Assuming no new shank if not at least 100 micron apart
ShankOpt = unique(Shank);
ShankID = nan(size(Shank));
for id = 1:length(ShankOpt)
    ShankID(Shank==ShankOpt(id))=id;
end
clusinfo.ShankID = ShankID;
Good_IDx = find(Good_ID);
