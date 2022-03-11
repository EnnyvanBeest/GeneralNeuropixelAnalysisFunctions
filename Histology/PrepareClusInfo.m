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
    Label = char(length(cluster_id));
    KSLabelfile = tdfread(fullfile(myKsDir,'cluster_KSLabel.tsv'));
    Label(ismember(cluster_id,KSLabelfile.cluster_id)) = KSLabelfile.KSLabel(ismember(KSLabelfile.cluster_id,cluster_id));
    Good_ID = ismember(cellstr(Label),'good'); %Identify good clusters
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
    
    % make sure cluster_id's match
    myOtherClusFile = dir(fullfile(myKsDir,'cluster_group.tsv'));% 
    cluster_groupOri = tdfread(fullfile(myOtherClusFile(1).folder,myOtherClusFile(1).name));
  
    KSLabel = clusinfo.KSLabel;
    Label = clusinfo.group; % You want the group, not the KSLABEL!
    depth = clusinfo.depth;
    channel = clusinfo.ch;
    Good_ID = ismember(cellstr(Label),'good'); %Identify good clusters
    
    % Make sure these two match
    Label(~(ismember(cluster_id,cluster_groupOri.cluster_id)),:)
    cluster_groupOri.group(~(ismember(cluster_groupOri.cluster_id,cluster_id)),:)

    if ~((sum(ismember(cluster_groupOri.cluster_id,cluster_id))/length(cluster_groupOri.cluster_id))==1) && ~((sum(ismember(cluster_id,cluster_groupOri.cluster_id))/length(cluster_id))==1)
        disp('Not the same amount of clusters!')
        keyboard
    end
    
    
%     [~,sortidx] =sort(cluster_id);
%     [~,sortidx2] = sort(cluster_groupOri.cluster_id);
%     A = cellstr(clusinfo.group(sortidx,:));
%     B = cellstr(cluster_groupOri.group(sortidx2,:));
%     A(find(~cell2mat(arrayfun(@(X) strcmp(A{X},B{X}),1:length(sortidx),'UniformOutput',0))))
%     B(find(~cell2mat(arrayfun(@(X) strcmp(A{X},B{X}),1:length(sortidx),'UniformOutput',0))))

end

% Find shank
myClusFile = dir(fullfile(myKsDir,'channel_map.npy'));
channelmap = readNPY(fullfile(myClusFile(1).folder,myClusFile(1).name));

myClusFile = dir(fullfile(myKsDir,'channel_positions.npy'));
channelpos = readNPY(fullfile(myClusFile(1).folder,myClusFile(1).name));

% Shank options
xpos = channelpos(:,1);
Shank = round(xpos(channel+1)/200); % Assuming no new shank if not at least 200 micron apart
ShankOpt = unique(Shank);
ShankID = nan(size(Shank));
for id = 1:length(ShankOpt)
    ShankID(Shank==ShankOpt(id))=id;
end
clusinfo.ShankID = ShankID;
Good_IDx = find(Good_ID);
