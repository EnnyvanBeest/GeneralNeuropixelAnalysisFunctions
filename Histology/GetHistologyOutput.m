histodone=0;
histoflag = 0;
if ~exist(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe))
    mkdir(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe))
end

if exist(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,'HistoEphysAlignment.mat'))
    tmpfile = load(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,'HistoEphysAlignment.mat'));
    try
        Depth2AreaPerChannel = tmpfile.Depth2AreaPerChannel;
        Depth2AreaPerUnit = tmpfile.Depth2AreaPerUnit;
        histodone=1;
        disp('Data already aligned to histology')
        histoflag=1;
    catch ME
        histodone=0;
    end
end
if ~histodone %Just in case it's not done yet
    histofile = dir(fullfile(myKsDir,'channel_locations.json'));
    if isempty(histofile)
        histofile = dir(fullfile(myKsDir,'*.csv'));
        if length(histofile)>1
            keyboard
        end
        if isempty(histofile)
            
            histofile = dir(fullfile(myKsDir,'probe_ccf.mat'));
            
            
            if isempty(histofile)
                disp('No channel locations known yet, do histology first')
                histoflag=0;
            else
                
                % alignment with petersen probe output
                histoflag = 1;
                histinfo =load(fullfile(histofile(1).folder,histofile(1).name));
                fullfile(fullfile(histofile(1).folder,histofile(1).name))

                % Align ephys data with probe
                [Depth2AreaPerChannel, Depth2AreaPerUnit]  = alignatlasdata(histinfo,AllenCCFPath,sp,clusinfo,1,1,2);
                
            end
        else % Automatic alignment with brain globe output
            histoflag = 1;
            histinfo =readtable(fullfile(histofile(1).folder,histofile(1).name),'ReadVariableNames',1,'Delimiter',',');
            
            % If available; find actual depth in brain
            % corresponding to this:
            trackcoordinates = readNPY(fullfile(histofile(1).folder,strrep(histofile(1).name,'.csv','.npy')));
            fullfile(fullfile(histofile(1).folder,histofile(1).name))
            % Align ephys data with probe
            [Depth2AreaPerChannel, Depth2AreaPerUnit]  = alignatlasdata(histinfo,AllenCCFPath,sp,clusinfo,1,0,2,trackcoordinates);
            
            
        end
    else
        histoflag=1;
        histinfo = fileread(fullfile(histofile(1).folder,histofile(1).name)); %Read json file
        histinfo = jsondecode(histinfo);% Decode json text
        fullfile(fullfile(histofile(1).folder,histofile(1).name))

        
        % Align ephys data with probe
        [Depth2AreaPerChannel, Depth2AreaPerUnit]  = alignatlasdata(histinfo,AllenCCFPath,sp,clusinfo,1,0,1);
    end
end
if histoflag
saveas(gcf,fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,'HistoEphysAlignment.fig'))
save(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,'HistoEphysAlignment.mat'),'Depth2AreaPerChannel','Depth2AreaPerUnit')
end