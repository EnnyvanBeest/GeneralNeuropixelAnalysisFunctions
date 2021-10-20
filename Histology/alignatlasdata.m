function [Depth2AreaPerChannel, Depth2AreaPerUnit] = alignatlasdata(histinfo,AllenCCFPath,sp,clusinfo,goodonly,surfacefirst,treeversion,trackcoordinates)
% Enny van Beest, based on AP_histology from AJPeters & IBLAPP from Mayo
% Faulkner
%% Inputs:
% histology info: probe track coordinates. Either output from
% AP_Histology pipeline (probe_ccf) or Brain Globe output (CSV file), 100a
% equally spaced points along probe track.
% AllenCCFPath: Path to AllenCCF (Github repository)
% Output from sp = loadKSdir(myKsDir); (Nick Steinmetz: https://github.com/cortex-lab/spikes)
% cluster information (KS/Phy output): channel (ID per cluster) and depth (also per Cluster)
%% Optional inputs: 
% surfacefirst = 1: position with lowest index is the surface of the brain. Default zero: Position with heighest index deeper in the brain
% treeversion: which Allen Brain Tree structure to use? (default 2, = 2017; 1 = older)
% trackcoordinates: for as many datapoints as in histinfo the X/Y/Z
% coordinates of the probe (e.g. the npy file from Brain Globe Output,
% using readNPY(fullfile(histofile(1).folder,strrep(histofile(1).name,'.csv','.npy')))



%% Use templatePositionsAmplitudes from the Spikes toolbox
[spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
    templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);
spikeCluster = sp.clu;
spikeTimes = sp.st;

%% Extract cluster info
try
    cluster_id = clusinfo.id;
catch
    cluster_id = clusinfo.cluster_id;
end
KSLabel = clusinfo.KSLabel;
depth = clusinfo.depth;
channel = clusinfo.ch;
Good_ID = find(sum(ismember(KSLabel,'good'),2)==4);

%% Select only good units to clean up MUA
if nargin>4 && goodonly
    spikeID = ismember(spikeCluster,Good_ID);
else
    spikeID = true(length(spikeCluster),1);
end

%% Surface first?
if nargin<6
    surfacefirst = 0;
end

%% Load structure tree allen brain
if nargin<7
    treeversion = 2;
end
if treeversion == 1
    tmp = readtable(fullfile(AllenCCFPath,'structure_tree_safe.csv'));
elseif treeversion == 2
    tmp = readtable(fullfile(AllenCCFPath,'structure_tree_safe_2017.csv'));
end
acronyms = lower(tmp.acronym);
color_hex = tmp.color_hex_triplet;
structure_id_path = tmp.structure_id_path;
id = tmp.id;

%% Actual coordinates known? - coordinates of track in Allen Brain space
coordinateflag =0;
if nargin>7
    coordinateflag = 1;
    X_ave=mean(trackcoordinates,1);            % mean; line of best fit will pass through this point
    dX=bsxfun(@minus,trackcoordinates,X_ave);  % residuals
    C=(dX'*dX)/(size(trackcoordinates,1)-1);           % variance-covariance matrix of X
    [R,D]=svd(C,0);             % singular value decomposition of C; C=R*D*R'
    
    D=diag(D);
    R2=D(1)/sum(D);
    disp(['Linear fit R2 = ' num2str(round(R2*1000)/10) '%'])
    
    x=dX*R(:,1);    % project residuals on R(:,1)
    x_min=min(x);
    x_max=max(x);
    dx=x_max-x_min;
    Xa=(x_min-0.05*dx)*R(:,1)' + X_ave;
    Xb=(x_max+0.05*dx)*R(:,1)' + X_ave;
    X_end=[Xa;Xb];
    
    figure
    plot3(X_end(:,1),X_end(:,2),X_end(:,3),'-r','LineWidth',3) % best fit line
    hold on
    plot3(trackcoordinates(:,1),trackcoordinates(:,2),trackcoordinates(:,3),'.k','MarkerSize',13)           % simulated noisy data
end

%% Open gui figure
gui_fig = figure('color','w');
flag =0; %to keep track of finishing this loop

%     if exist('gui_fig')
%         delete(gui_fig)
%     end
% subplot(3,5,[1:2,6:7,11:12])
% plotDriftmap(spikeTimesCorrected, spikeAmps, spikeDepths);
% title('DriftMap')
depthunique = unique(spikeDepths(spikeID));

% Calculate number of spikes across depth
nrspikesperdepth = cell2mat(arrayfun(@(X) sum(spikeTimes(spikeDepths==X & spikeID==1)),depthunique,'UniformOutput',0));
thresh = quantile(nrspikesperdepth,0.01);

%Find 'gaps' of low activity:
gaps = depthunique(find(nrspikesperdepth<thresh));
endpoint = max(depthunique);
startpoint = min(depthunique);

subplot(3,5,[1,6,11])
% Calculate number of spikes across depth
avgAmplitude = cell2mat(arrayfun(@(X) nanmean(spikeAmps(spikeDepths==X & spikeID==1)),depthunique,'UniformOutput',0));
box off
plot(smooth(avgAmplitude),depthunique,'-')
title('avg Amplitude')
ylim([startpoint,endpoint]);

box off
subplot(3,5,[2,7,12])
plot(smooth(nrspikesperdepth),depthunique,'-'); hold on;
title('Nr Spikes')
line([thresh thresh],get(gca,'ylim'),'color',[1 0 0],'LineWidth',2)
ylim([startpoint,endpoint]);
hold off

%% Get multiunit correlation - Copied from Petersen github
n_corr_groups = 40;
depth_group_edges = linspace(startpoint,endpoint,n_corr_groups+1);
depth_group = discretize(spikeDepths,depth_group_edges);
depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);
unique_depths = 1:length(depth_group_edges)-1;

spike_binning = 0.01; % seconds
corr_edges = nanmin(spikeTimes(spikeID==1)):spike_binning:nanmax(spikeTimes(spikeID==1));
corr_centers = corr_edges(1:end-1) + diff(corr_edges);

binned_spikes_depth = zeros(length(unique_depths),length(corr_edges)-1);
for curr_depth = 1:length(unique_depths)
    binned_spikes_depth(curr_depth,:) = histcounts(spikeTimes(depth_group == unique_depths(curr_depth) & spikeID==1), corr_edges);
end

mua_corr = corrcoef(binned_spikes_depth');
limup = quantile(mua_corr(:),0.8);
% Plot multiunit correlation
multiunit_ax = subplot(3,5,[3:4,8:9,13:14]);
h=imagesc(depth_group_centers,depth_group_centers,mua_corr);
caxis([0,max(mua_corr(mua_corr~=1))]); colormap(hot);
set(h,'Alphadata',~isnan(mua_corr))
set(gca,'Color',[0.5 0.5 0.5])
ylim([startpoint,endpoint]);
xlim([startpoint,endpoint]);
set(multiunit_ax,'YDir','normal');
title('MUA correlation');
xlabel(multiunit_ax,'Multiunit depth');

while ~flag
    %Now divide position of probe along this track
    if istable(histinfo)&& any(histinfo.Position)
        histinfo.RegionAcronym(ismember(histinfo.RegionAcronym,'Not found in brain')) = {'root'};
        if ~surfacefirst
            areapoints = linspace(startpoint,endpoint,length(histinfo.Position));
            if coordinateflag
                trackcoordinates = [linspace(X_end(1,1),X_end(2,1),length(areapoints));linspace(X_end(1,2),X_end(2,2),length(areapoints));linspace(X_end(1,3),X_end(2,3),length(areapoints))]';
            end
        else
            areapoints = linspace(endpoint,startpoint,length(histinfo.Position));
            if coordinateflag
                trackcoordinates = [linspace(X_end(2,1),X_end(1,1),length(areapoints));linspace(X_end(2,2),X_end(1,2),length(areapoints));linspace(X_end(2,3),X_end(1,3),length(areapoints))]';
            end
        end
        [UniqueAreas,IA,IC] = unique((histinfo.RegionAcronym),'stable');
    elseif isstruct(histinfo)&& isfield(histinfo,'probe_ccf')
        if ~surfacefirst
            areapoints = (linspace(startpoint,endpoint,length(histinfo.probe_ccf.trajectory_coords)));
        else
            areapoints = (linspace(endpoint,startpoint,length(histinfo.probe_ccf.trajectory_coords)));
        end
        histinfo.probe_ccf.trajectory_acronyms = acronyms(histinfo.probe_ccf.trajectory_areas);
        [UniqueAreas,IA,IC] = unique(histinfo.probe_ccf.trajectory_acronyms,'stable');
        histinfo.RegionAcronym = histinfo.probe_ccf.trajectory_acronyms;
    else
        areapoints = nan(1,max(channel)+1);
        histinfo.RegionAcronym  = cell(1,max(channel)+1);
        for chid = 1:max(channel)+1
            eval(['areapoints(1,' num2str(chid) ')= abs(histinfo.channel_' num2str(chid-1) '.z);'])
            eval([' histinfo.RegionAcronym {1,' num2str(chid) '}= histinfo.channel_' num2str(chid-1) '.brain_region;'])
        end
        % not always in the correct order, align
        if sum(ismember([-1,1],unique(sign(diff(areapoints)))))==2
            [areapoints, sortid] = sort(areapoints,'descend');
            histinfo.RegionAcronym = histinfo.RegionAcronym(sortid);
        end
        histinfo.RegionAcronym(ismember(histinfo.RegionAcronym,'void')) = {'root'};
        
        [UniqueAreas,IA,IC] = unique(fliplr(histinfo.RegionAcronym),'stable');
      
            
        
    end
    UniqueAreas = lower(UniqueAreas); % case insensitive
    switchpoints = [1; find(diff(IC)~=0)+1; length(IC)]; %Find points where region changes
    AllAreas = histinfo.RegionAcronym(switchpoints);
    
    if exist('probe_areas_ax')
        delete(probe_areas_ax)
    end
    probe_areas_ax  = subplot(3,5,[5,10,15]);
    patchobj = gobjects;
    textobj = gobjects;
    for i=2:length(switchpoints)
        patchobj(i-1) = patch([0 1 1 0],[areapoints(switchpoints(i-1)) areapoints(switchpoints(i-1)) areapoints(switchpoints(i)) areapoints(switchpoints(i))],hex2rgb(color_hex(ismember(acronyms,UniqueAreas{IC(switchpoints(i-1))}))));
        textobj(i-1) = text(0.5,nanmean([areapoints(switchpoints(i-1)) areapoints(switchpoints(i))]),UniqueAreas{IC(switchpoints(i-1))},'HorizontalAlignment','center');
    end
    patchobj(i) = patch([0 1 1 0],[areapoints(switchpoints(i)) areapoints(switchpoints(i)) areapoints(end) areapoints(end)],hex2rgb(color_hex(ismember(acronyms,UniqueAreas{IC(switchpoints(i))}))));
    textobj(i) = text(0.5,nanmean([areapoints(switchpoints(i)) areapoints(end)]),UniqueAreas{IC(switchpoints(i))},'HorizontalAlignment','center');
    ylim([startpoint,endpoint]);
    title({'Probe areas','(ws keys to move, a to add ref line, d to delete ref line, f for flip probe ori, 123 for factor)','(q: save & quit, r: reset)'});
    if coordinateflag
        
        yyaxis right
        ylim([startpoint,endpoint]);
        %Find corresponding trackcoordinates
        tmplabel=trackcoordinates(cell2mat(arrayfun(@(X) find(abs(areapoints-X)==min(abs(areapoints-X)),1,'first'),get(gca,'YTick'),'UniformOutput',0)),:);
        tmplabel = num2cell(tmplabel,2);
        tmplabel = cellfun(@(X) ['[' num2str(round(X(1))) ';', num2str(round(X(2))), ';', num2str(round(X(3))),']'],tmplabel,'UniformOutput',0);
        set(gca,'YTickLabel',tmplabel)
        yyaxis left
    end
    
    % Draw corresponding area lines on multi unit
    subplot(multiunit_ax)
    % Draw boundary lines at borders (and undo clipping to extend across all)
    if exist('boundary_lines')
        delete(boundary_lines)
    end
    boundary_lines = gobjects;
    for curr_boundary = 1:length(switchpoints)
        boundary_lines(curr_boundary,1) = line(probe_areas_ax,[0 1], ...
            repmat(areapoints(switchpoints(curr_boundary)),1,2),'color','b','linewidth',1);
        boundary_lines(curr_boundary,2) = line(multiunit_ax,[startpoint endpoint], ...
            repmat(areapoints(switchpoints(curr_boundary)),1,2),'color','y','linewidth',1,'LineStyle','--');
    end
    %% Interface
    matchedswitchpoints = nan(2,length(switchpoints));
    matchedswitchpoints(1,:)=areapoints(switchpoints);
    newswitchpoints = switchpoints;
    if coordinateflag
        newtrackcoordinates = nan(size(trackcoordinates));
    end
    newareapoints = areapoints; %new s
    oristartpoint = startpoint;
    oriendpoint = endpoint;
    stepsize = unique(round(abs(diff(areapoints))));
    disp('ws keys to move, a to add reference line, d to delete ref line, f for flip probe orientation, 123 for speed of movement, q: save & quit, r: reset');
    okay = 0;
    key = '';
    y_change = 3;
    while ~okay
        switch key
            % Set amounts to move by with/without shift
            case '1'
                y_change = 1;
            case '2'
                y_change = 10;
            case '3'
                y_change = 100;
                % up/down: move probe areas
            case 'w'
                if isnan( matchedswitchpoints(2,[1 size(matchedswitchpoints,2)]))
                    matchedswitchpoints(2,[1 size(matchedswitchpoints,2)]) = matchedswitchpoints(1,[1 size(matchedswitchpoints,2)]) - y_change;
                else
                    matchedswitchpoints(2,[1 size(matchedswitchpoints,2)]) = matchedswitchpoints(2,[1 size(matchedswitchpoints,2)]) - y_change;
                end
            case 's'
                if isnan(matchedswitchpoints(2,[1 size(matchedswitchpoints,2)]))
                    matchedswitchpoints(2,[1 size(matchedswitchpoints,2)]) = matchedswitchpoints(1,[1 size(matchedswitchpoints,2)]) + y_change;
                else
                    matchedswitchpoints(2,[1 size(matchedswitchpoints,2)]) = matchedswitchpoints(2,[1 size(matchedswitchpoints,2)]) + y_change;
                end
            case 'r'
                newswitchpoints = switchpoints;
                matchedswitchpoints(2,:) = nan;
                for curr_boundary = 1:length(newswitchpoints)
                    set(boundary_lines(curr_boundary,1),'color','b')
                    set(boundary_lines(curr_boundary,2),'color','y')
                end
            case 'a' %Add reference line
                disp('Click to add reference line on probe')
                roi1 = drawpoint;
                [~,minidx] = nanmin(abs(areapoints(newswitchpoints)-roi1.Position(2)));
                set(boundary_lines(minidx,1),'color','r')
                delete(roi1)
                %find closest line;
                disp('Click to add reference line on MUA')
                roi2 = drawpoint;
                set(boundary_lines(minidx,2),'color','r','YData',[roi2.Position(2) roi2.Position(2)])
                matchedswitchpoints(2,minidx) = roi2.Position(2);
                delete(roi2)
            case 'd' %Remove reference line
                disp('Click to remove reference line on probe')
                roi1 = drawpoint;
                [~,minidx] = nanmin(abs(areapoints(newswitchpoints)-roi1.Position(2)));
                set(boundary_lines(minidx,1),'color','b')
                delete(roi1)
                set(boundary_lines(minidx,2),'color','y','YData',[areapoints(newswitchpoints(minidx)) areapoints(newswitchpoints(minidx))])
                matchedswitchpoints(2,minidx) = nan;
                delete(roi2)
                % q: save and quit
            case 'f' %Flip orientation of probe
                surfacefirst = abs(surfacefirst-1);
                break
            case 'q'
                okay = 1;
                flag = 1;
                break
        end
        newswitchpoints(newswitchpoints<1) = nan;
        newswitchpoints(newswitchpoints>length(newareapoints))=nan;
        % Update figure
        if sum(~isnan(matchedswitchpoints(2,:))) > 1
            %         match the two
            nonnanidx = find(~isnan(matchedswitchpoints(2,:)));
            if nonnanidx(1)~=1
                %         1 to first matchedswitchpoint
                newvals = matchedswitchpoints(2,1:nonnanidx(1));
                oldvals = matchedswitchpoints(1,1:nonnanidx(1));
                proportion = (oldvals-oldvals(1))./(oldvals(end)-oldvals(1)); %proportion of areas in between
                %New switchpoints, keep proportion of in between areas the same
                newswitchpoints(1:nonnanidx(1)) = cell2mat(arrayfun(@(X) find(abs(newareapoints-X)==nanmin(abs(newareapoints-X)),1,'first'),(newvals(end)-oldvals(1))*proportion+oldvals(1),'UniformOutput',0));
                if coordinateflag
                    newtrackcoordinates(newswitchpoints(1):newswitchpoints(nonnanidx(1)),:)= [linspace(trackcoordinates(switchpoints(1),1),trackcoordinates(switchpoints(nonnanidx(1)),1),length(newswitchpoints(1):newswitchpoints(nonnanidx(1)))); ...
                        linspace(trackcoordinates(switchpoints(1),2),trackcoordinates(switchpoints(nonnanidx(1)),2),length(newswitchpoints(1):newswitchpoints(nonnanidx(1)))); ...
                        linspace(trackcoordinates(switchpoints(1),3),trackcoordinates(switchpoints(nonnanidx(1)),3),length(newswitchpoints(1):newswitchpoints(nonnanidx(1))))]';
                end
            end
            for i = 1:length(nonnanidx)-1
                newvals = matchedswitchpoints(2,nonnanidx(i):nonnanidx(i+1));
                oldvals = matchedswitchpoints(1,nonnanidx(i):nonnanidx(i+1));
                proportion = (oldvals-oldvals(1))./(oldvals(end)-oldvals(1)); %proportion of areas in between
                %New switchpoints, keep proportion of in between areas the same
                newswitchpoints(nonnanidx(i):nonnanidx(i+1)) = cell2mat(arrayfun(@(X) find(abs(newareapoints-X)==nanmin(abs(newareapoints-X)),1,'first'),(newvals(end)-newvals(1))*proportion+newvals(1),'UniformOutput',0));
                if coordinateflag
                    newtrackcoordinates(newswitchpoints(nonnanidx(i)):newswitchpoints(nonnanidx(i+1)),:)= [linspace(trackcoordinates(switchpoints(nonnanidx(i)),1),trackcoordinates(switchpoints(nonnanidx(i+1)),1),length(newswitchpoints(nonnanidx(i)):newswitchpoints(nonnanidx(i+1)))); ...
                        linspace(trackcoordinates(switchpoints(nonnanidx(i)),2),trackcoordinates(switchpoints(nonnanidx(i+1)),2),length(newswitchpoints(nonnanidx(i)):newswitchpoints(nonnanidx(i+1))));...
                        linspace(trackcoordinates(switchpoints(nonnanidx(i)),3),trackcoordinates(switchpoints(nonnanidx(i+1)),3),length(newswitchpoints(nonnanidx(i)):newswitchpoints(nonnanidx(i+1))))]';
                end
            end
            %Now the bit after
            if nonnanidx(i+1)<size(matchedswitchpoints,2)
                newvals = matchedswitchpoints(2,nonnanidx(i+1):end);                
                oldvals = matchedswitchpoints(1,nonnanidx(i+1):end);
                proportion = ((oldvals-oldvals(1))./(oldvals(end)-oldvals(1))); %proportion of areas in between
                %New switchpoints, keep proportion of in between areas the same
                newswitchpoints(nonnanidx(i+1):end) = cell2mat(arrayfun(@(X) find(abs(newareapoints-X)==nanmin(abs(newareapoints-X)),1,'first'),(oldvals(end)-newvals(1))*proportion+newvals(1),'UniformOutput',0));
                if coordinateflag
                    newtrackcoordinates(newswitchpoints(nonnanidx(i+1)):length(newtrackcoordinates),:)= [linspace(trackcoordinates(switchpoints(nonnanidx(i+1)),1),trackcoordinates(end,1),length(newswitchpoints(nonnanidx(i+1)):length(newtrackcoordinates))); ...
                        linspace(trackcoordinates(switchpoints(nonnanidx(i+1)),2),trackcoordinates(end,2),length(newswitchpoints(nonnanidx(i+1)):length(newtrackcoordinates))); ...
                        linspace(trackcoordinates(switchpoints(nonnanidx(i+1)),3),trackcoordinates(end,3),length(newswitchpoints(nonnanidx(i+1)):length(newtrackcoordinates)))]';
                end
            end
            
            if any(diff(newswitchpoints)<0)
                disp('Something went wrong, Reset')
                newswitchpoints = switchpoints;
                matchedswitchpoints(2,:) = nan;
                for curr_boundary = 1:length(newswitchpoints)
                    set(boundary_lines(curr_boundary,1),'color','b')
                    set(boundary_lines(curr_boundary,2),'color','y')
                end
            end
            for i=2:length(newswitchpoints)
                patchobj(i-1).YData=[newareapoints(newswitchpoints(i-1)) newareapoints(newswitchpoints(i-1)) newareapoints(newswitchpoints(i)) newareapoints(newswitchpoints(i))];
                textobj(i-1).Position(2) = nanmean([newareapoints(newswitchpoints(i-1)) newareapoints(newswitchpoints(i))]);
            end
            patchobj(i).YData=[newareapoints(newswitchpoints(i)) newareapoints(newswitchpoints(i)) newareapoints(newswitchpoints(end)) newareapoints(newswitchpoints(end))];
            textobj(i).Position(2) = nanmean([newareapoints(newswitchpoints(i)) newareapoints(newswitchpoints(end))]);
            
            ylim([oristartpoint oriendpoint])
            
            % update boundary lines at borders (and undo clipping to extend across all)
            for curr_boundary = 1:length(newswitchpoints)
                set(boundary_lines(curr_boundary,1),'YData',repmat(newareapoints(newswitchpoints(curr_boundary)),1,2))
                set(boundary_lines(curr_boundary,2),'YData',repmat(newareapoints(newswitchpoints(curr_boundary)),1,2))
            end
            
            if coordinateflag
                subplot(probe_areas_ax)
                yyaxis right
                hold on
                
                ylim([oristartpoint,oriendpoint]);
                %Find corresponding trackcoordinates
                tmplabel=newtrackcoordinates(cell2mat(arrayfun(@(X) find(abs(newareapoints-X)==min(abs(newareapoints-X)),1,'first'),get(gca,'YTick'),'UniformOutput',0)),:);
                tmplabel = num2cell(tmplabel,2);
                tmplabel = cellfun(@(X) ['[' num2str(round(X(1))) ';', num2str(round(X(2))), ';', num2str(round(X(3))),']'],tmplabel,'UniformOutput',0);
                set(gca,'YTickLabel',tmplabel)
                yyaxis left
                
            end
        end
        
        %Input?
        waitforbuttonpress
        key = get(gcf,'CurrentCharacter');
    end
end
%% Shift area
areasaligned = cell(1,length(histinfo.RegionAcronym));
for i = 1:length(newswitchpoints)-1
    areasaligned(newswitchpoints(i):newswitchpoints(i+1)) = lower(AllAreas(i));
end
areasaligned(cell2mat(cellfun(@isempty,areasaligned,'UniformOutput',0)))={'root'};

%% Make a depth to area conversion table - per channel
[channels,IA,IC]=unique(channel);
% depth = depth(IA);
if coordinateflag
    Depth2AreaPerChannel=cell(4,max(channels)+1); %Depth, regionname, color, actual coordinates
else
    Depth2AreaPerChannel=cell(3,max(channels)+1); %Depth, regionname, color
end
Depth2AreaPerChannel(1,channels+1) = (arrayfun(@(X) nanmean(depth(IC==X)),1:length(channels),'UniformOutput',0)); % on the probe
tmp = (cellfun(@(X) find(abs(areapoints-X)==min(abs(areapoints-X)),1,'first'),(arrayfun(@(X) nanmean(depth(IC==X)),1:length(channels),'UniformOutput',0)),'UniformOutput',0));
idx = channels+1;
idx(cell2mat(cellfun(@isempty,tmp,'UniformOutput',0)))=[];
Depth2AreaPerChannel(2,idx) = areasaligned(cell2mat(tmp));
Depth2AreaPerChannel(3,idx) = cellfun(@(X) color_hex(ismember(acronyms,X)),Depth2AreaPerChannel(2,idx),'UniformOutput',0);
if coordinateflag
    Depth2AreaPerChannel(4,idx) = num2cell(newtrackcoordinates(cell2mat(tmp),:),2);
end
%% Make a depth to area conversion table - per unit
tmp = (arrayfun(@(X) find(abs(areapoints-X)==min(abs(areapoints-X)),1,'first'),depth,'UniformOutput',0));
clusterarea = repmat({'unknown'},1,length(cluster_id));
idx = 1:length(clusterarea);
idx(cell2mat(cellfun(@isempty,tmp,'UniformOutput',0))) = [];
clusterarea(idx) = areasaligned(cell2mat(tmp(idx)));
clustercolor = repmat({'#808080'},1,length(cluster_id));
clustercolor(idx) = cellfun(@(X) color_hex(ismember(acronyms,X)),clusterarea(idx),'UniformOutput',0);

if coordinateflag
    clustercoord = repmat({[nan,nan,nan]},1,length(cluster_id));
    clustercoord(idx) = num2cell(newtrackcoordinates(cell2mat(tmp(idx)),:),2);
    Depth2AreaPerUnit = table(cluster_id,depth,clusterarea',clustercolor',clustercoord','VariableNames',{'Cluster_ID','Depth','Area','Color','Coordinates'});
else
    Depth2AreaPerUnit = table(cluster_id,depth,clusterarea',clustercolor','VariableNames',{'Cluster_ID','Depth','Area','Color'});
end

end







