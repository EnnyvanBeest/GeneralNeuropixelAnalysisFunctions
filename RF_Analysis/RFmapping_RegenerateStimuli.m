%% Load and save RF mapping

% For mpep experiments;
DataDir = '\\128.40.198.18\Subjects\';%'Z:\' %
rigname = 'Zatteo';
MiceOpt = {'EB010','EB011'};%{'FT039','CB020','Charu_Inscopix1','EB001','EB002','EB003','EB004','EB005','EB006','EB007','EB008','EB009','CB007','CB008'}; %,'EB001','EB003'}%{'CB008'};%{'EB001'}%{'EB001','EB002','EB003','CB007','CB008'};%,'CB007','CB008'} %'CB007'

DateOpt = cellfun(@(X) dir(fullfile(DataDir,X,'*-*')),MiceOpt,'UniformOutput',0);
DateOpt = cellfun(@(X) X([X.isdir]),DateOpt,'UniformOutput',0);
DateOpt = cellfun(@(X) {X.name},DateOpt,'UniformOutput',0);

%Rig info
RigInfo = RigInfoGet(rigname);
myScreenInfo = ScreenInfo(RigInfo);
myScreenInfo = myScreenInfo.CalibrationLoad;
myScreenInfo.windowPtr = nan;

%Loop over mice/ sessions
for midx = 1:length(MiceOpt)
    MiceOpt{midx}
    Dates4Mouse = DateOpt{midx};
    for didx = 1:length(Dates4Mouse)
        % Within folders, look for 'RF mapping sessions'
        thisdate = Dates4Mouse{didx}
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
                    continue
                end
                fileID = fopen(fullfile(listfiles(idx).folder,listfiles(idx).name));
                A = fscanf(fileID,'%c');
                fclose(fileID)
                
                if any(strfind(A,'SparseNoise'))
                    % Load protocol to check for contrast >0
                    try
                    Protocol = load(fullfile(fullfile(listfiles(idx).folder,'Protocol.mat')));
                    Protocol = Protocol.Protocol;
                    if Protocol.pars(strcmp(Protocol.parnames,'c')) == 0
                        %not a RF session, continue
                        continue
                    end
                    
                    RFsess = [RFsess {subsess(sesidx).name}];
                    flag = 1;
                    continue
                    catch 
                        continue
                    end
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
            if exist(fullfile(DataDir,MiceOpt{midx},thisdate,thisses,'RF_MappingStimuli.mat'))
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
            ntrials = length(protocol.seqnums);
            SS = cell(1,ntrials);
            flagseed=0;
            for trialid=1:ntrials
                if isnan(protocol.pars(ismember(protocol.parnames,'seed'))) || flagseed
                    %replace seed = nan with seq. number
                    flagseed=1;
                    protocol.pars(ismember(protocol.parnames,'seed')) = protocol.seqnums(trialid);                    
                end
                matlabfunctused = strrep(protocol.xfile,'.x','');
                
                %Run Script              
                tmp = feval(matlabfunctused,myScreenInfo,protocol.pars);
                
                %Screenstim cannot be saved. Save individual parameters
                %seperately
                SS{trialid} = struct;
                SS{trialid}.ImageTextures = tmp.ImageTextures;
                SS{trialid}.Parameters = tmp.Parameters;
                SS{trialid}.ImageSequence = tmp.ImageSequence;


            end
            if ~isempty(SS)
                save(fullfile(DataDir,MiceOpt{midx},thisdate,thisses,'RF_MappingStimuli.mat'),'SS')
                clear StimInformation
            end
        end
    end
end