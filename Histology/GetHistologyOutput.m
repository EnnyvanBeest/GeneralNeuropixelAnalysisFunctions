histodone=0;
histoflag = 0;
if ~exist(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe))
    mkdir(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe))
end

if exist(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,'HistoEphysAlignment.mat')) && ~NewHistologyNeeded
    tmpfile = load(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,'HistoEphysAlignment.mat'));
    try
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
            disp('Detecting multiple files... npix2 probe?')
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
                Depth2AreaPerUnit  = alignatlasdata(histinfo,AllenCCFPath,sp,clusinfo,1,1,2);
                
            end
        else % Automatic alignment with brain globe output
            histoflag = 1;
            histinfo = arrayfun(@(X) readtable(fullfile(histofile(X).folder,histofile(X).name),'ReadVariableNames',1,'Delimiter',','),1:length(histofile),'UniformOutput',0);
            
            % If available; find actual depth in brain
            % corresponding to this:
            trackcoordinates = arrayfun(@(X) readNPY(fullfile(histofile(X).folder,strrep(histofile(X).name,'.csv','.npy'))),1:length(histofile),'UniformOutput',0);
            % Align ephys data with probe
            Depth2AreaPerUnit  = alignatlasdata(histinfo,AllenCCFPath,sp,clusinfo,0,1,2,trackcoordinates);   
        end
    else
        histoflag=1;
        histinfo = fileread(fullfile(histofile(1).folder,histofile(1).name)); %Read json file
        histinfo = jsondecode(histinfo);% Decode json text
        fullfile(fullfile(histofile(1).folder,histofile(1).name))

        % Align ephys data with probe
        Depth2AreaPerUnit  = alignatlasdata(histinfo,AllenCCFPath,sp,clusinfo,1,0,1);
    end
end
if histoflag
saveas(gcf,fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,'HistoEphysAlignment.fig'))
save(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,'HistoEphysAlignment.mat'),'Depth2AreaPerUnit')
end