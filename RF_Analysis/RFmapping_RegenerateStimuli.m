%% Load and save RF mapping

% For mpep experiments;
DataDir = '\\znas\Subjects'
rigname = 'Zatteo';
MiceOpt = {'EB001','EB002','EB003','EB004','EB005'}%,'CB007','CB008'}; %,'EB001','EB003'}%{'CB008'};%{'EB001'}%{'EB001','EB002','EB003','CB007','CB008'};%,'CB007','CB008'} %'CB007'

DateOpt = cellfun(@(X) dir(fullfile(DataDir,X,'*-*')),MiceOpt,'UniformOutput',0);
DateOpt = cellfun(@(X) X([X.isdir]),DateOpt,'UniformOutput',0);
DateOpt = cellfun(@(X) {X.name},DateOpt,'UniformOutput',0);

for midx = 1:length(MiceOpt)
    Dates4Mouse = DateOpt{midx};
    for didx = 1:length(Dates4Mouse)
        % Within folders, look for 'RF mapping sessions'
        thisdate = Dates4Mouse{didx};
        subsess = dir(fullfile(DataDir,MiceOpt{midx},Dates4Mouse{didx}));
        subsess(1:2) = []; %remove '.' and '..'
        flag = 0;
        RFsess = [];
        for sesidx=1:length(subsess)
            listfiles = dir(fullfile(subsess(sesidx).folder,subsess(sesidx).name));
            listfiles(~cell2mat(cellfun(@isempty,cellfun(@(X) strfind(X,'.pickle'),{listfiles(:).name},'UniformOutput',0),'UniformOutput',0))) = [];
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
        for sesidx=1:length(RFsess)
            thisses=RFsess{sesidx}
            if 0%exist(fullfile(DataDir,MiceOpt{midx},thisdate,thisses,'RF_MappingStimuli.mat'))
                disp([fullfile(DataDir,MiceOpt{midx},thisdate,thisses,'RF_MappingStimuli.mat') ' already exists, skip...'])
                continue
            end
            protocoldfile = dir(fullfile(DataDir,MiceOpt{midx},thisdate,thisses,'Protocol.mat'));
            if isempty(protocoldfile)
                disp('No RF mapping file found... skip')
                continue
            end
            protocol = load(fullfile(protocoldfile.folder,protocoldfile.name));
            protocol = protocol.Protocol;
            matlabfunctused = strrep(protocol.xfile,'.x','');
            
            %Run Script
            RigInfo = RigInfoGet(rigname);
            myScreenInfo = ScreenInfo(RigInfo);
            myScreenInfo = myScreenInfo.CalibrationLoad;
            myScreenInfo.windowPtr = nan;
            SS = feval(matlabfunctused,myScreenInfo,protocol.pars);
            
            if ~isempty(SS)
                save(fullfile(DataDir,MiceOpt{midx},thisdate,thisses,'RF_MappingStimuli.mat'),'SS')
                clear SS
            end
        end
    end
end