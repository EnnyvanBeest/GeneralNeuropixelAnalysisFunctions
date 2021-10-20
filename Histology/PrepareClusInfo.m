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
    Good_ID = find(KSLabel=='g');
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
    save(fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,thisprobe,'CuratedResults.mat'),'CurationDone')
    clusinfo = tdfread(fullfile(myClusFile(1).folder,myClusFile(1).name));
    curratedflag=1;
    
    cluster_id = clusinfo.id;
    KSLabel = clusinfo.KSLabel;
    depth = clusinfo.depth;
    channel = clusinfo.ch;
    Good_ID = find(sum(ismember(KSLabel,'good'),2)==4);
    
end

Good_IDx = clusinfo.id(Good_ID);
