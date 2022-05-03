%% User Input
DataDir = {'\\znas\Subjects','\\128.40.198.18\Subjects'}%'\\znas\Subjects' %'
SaveDir = 'E:\Data\Results_Apr2022\'
tmpdatafolder = 'D:\tmpdata\';
LocalDir = 'E:\Data\PyKSOutput';% 'E:\Data\KiloSortOutput';%
AllenCCFPath = 'C:\Users\EnnyB\Documents\MATLAB\allenCCF'
HistoFolder = 'E:\Data\Histology\';
SaveFiguresTo = 'E:\Data\Figures'
MiceOpt = {'EB001','EB002','EB003','EB004','EB005','EB006','EB007','EB008','EB009','EB010','EB011','CB020'};%{'EB001','EB002','EB003','EB004','EB005','EB006','EB007','EB008','EB009'};% CB020;%{};%{}%{'EB006'};%,'CB007','CB008'}; %,'EB001','EB003'}%{'CB008'};%{'EB001'}%{'EB001','EB002','EB003','CB007','CB008'};%,'CB007','CB008'} %'CB007'
nidq_sync_used = zeros(1,length(MiceOpt));
nidq_sync_used(ismember(MiceOpt,{'EB001','CB007','CB008'}))=1;
DataDir2Use = repmat(2,[1,length(MiceOpt)]);
DataDir2Use(ismember(MiceOpt,{'EB001','EB002','EB003','EB004','EB005','CB007','CB008'}))=1;
ProbeType = repmat({'1_3b'},1,length(MiceOpt)); % if multiple probe types, specify which animals have 2_4s
ProbeType(ismember(MiceOpt,{'FT039','CB020','EB014'}))={'2_4S'};

%% Add all subdirectories
addpath(genpath(cd))

%% Get histology alignment
MainHistology

%% Run mpep analysis
Generalmpep_Analysis