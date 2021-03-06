%% This script is to find receptive fields after presentation of sparse noise (i.e. mpep 'asyncsparsenoise.m')
%% Before doing analysis, run 'RFmapping_RegenerateStimuli' on the computer where you ran the mpep asyncsparsenoise script to get the exact presentation pattersn


%% User Input (you probably won't need all of these, but see RF script)
DataDir = {'\\znas\Subjects','\\128.40.198.18\Subjects'}%'\\znas\Subjects' %'
SaveDir = 'E:\Data\ResultsOngoing\'
tmpdatafolder = 'D:\tmpdata\';
BackUpPath = 'F:\' %Only for back-up
LocalDir = 'E:\Data\PyKSOutput';% 'E:\Data\KiloSortOutput';%
AllenCCFPath = 'C:\Users\EnnyB\Documents\MATLAB\allenCCF'
storevideopath=fullfile(tmpdatafolder,'Videos');
HistoFolder = 'E:\Data\Histology\';
SaveFiguresTo = 'E:\Data\Figures'
MiceOpt = {'EB001','EB003','EB006','EB007','EB008','EB009','EB010','EB011','EB012','EB013'}%{'EB014'};%{'EB016','EB017','EB018','EB019','AV009'};%,'EB014'}% {};%{}; %{'EB009','EB010','EB011'};%{};%{'EB007','EB008','EB009'};%{'EB001','EB002','EB003','EB004','EB005','EB006','EB007','EB008','EB009','EB010','EB011','CB020'};%{'EB001','EB002','EB003','EB004','EB005','EB006','EB007','EB008','EB009'};% CB020;%{};%{}%{'EB006'};%,'CB007','CB008'}; %,'EB001','EB003'}%{'CB008'};%{'EB001'}%{'EB001','EB002','EB003','CB007','CB008'};%,'CB007','CB008'} %'CB007'
nidq_sync_used = zeros(1,length(MiceOpt));
nidq_sync_used(ismember(MiceOpt,{'EB001','CB007','CB008'}))=1;
DataDir2Use = repmat(2,[1,length(MiceOpt)]);
DataDir2Use(ismember(MiceOpt,{'EB001','EB002','EB003','EB004','EB005','CB007','CB008','AL056'}))=1;
ProbeType = repmat({'1_3b'},1,length(MiceOpt));
ProbeType(ismember(MiceOpt,{'FT039','CB020','EB014'}))={'2_4S'};
VideoTypes = {'bellyCam','eyeCamRight','eyeCamLeft','eye'};

RunQualityMatrix = 0;
InspectQualityMatrix = 1;
maxsessnr = 2; %max nr. sessions on a day (doesn't need to be accurate)
MinDist2Include = 85; %Including only trials in which mouse reached at least xcm for some statistics (past all stimuli?)
pretrialtime = 2; %take up to x seconds prior trial
posttrialtime = 2; % take up to x seconds post trial
addpath(genpath(cd))
addpath(genpath('C:\Users\EnnyB\Documents\GitHub\spikes'))
addpath(genpath('C:\Users\EnnyB\Documents\GitHub\npy-matlab'))
addpath(genpath('C:\Users\EnnyB\Documents\GitHub\AP_histology'))
addpath(genpath('C:\Users\EnnyB\Documents\GitHub\MainAnalysis_SpatialNavigation'))

%% RF mapping analysis
LoadSpikesRFBoot

