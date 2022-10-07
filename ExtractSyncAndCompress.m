%% Compress & Move Neuropixels data to server
% Enny van Beest, Oct 2022
% Loop over local bin files to:
% 1. Extract sync
% 2. Bin --> cbin
% 3. Upload everything (.cbin, meta and sync file) except raw .bin to server
% 4. Shut down computer (optionally)
% OR, in the case the .cbin file already exists locally
% 1. Check whether data is copied to server already (compare sizes etc.)
% 2. Check when this was done. if >xclone days, delete the local copy
% After looping, shut down the computer optionally
LocalDir = 'E:\SpikeGLX\' % This is your local directory. The assumption is that all data here will have to be compressed and moved to the server
ServerDir = '\\zaru.cortexlab.net\Subjects\'
XCloneDays = 3; % The number of days on the server after which we can assume a clone has been made
TurnOffDesktop = 1; % Turn off this PC after running the script

% Find all folders with bin files
localEphysFiles = dir(fullfile(LocalFolder,'**','*.ap.bin'));
if isempty(localEphysFiles)
    fprintf('There are no ephys files in the local directory...')
    return
end

for fid = 1:length(localEphysFiles)
    disp(['This is ' localEphysFiles(fid).name])
    % Extract full path on server from local Ephys File name
    fileParts = strsplit(localEphysFiles(fid).name,'_');
    % Find a date in this (assuming yyyy-mm-dd)
    thisDate = fileParts{find(cell2mat(cellfun(@(X) sum(ismember(X,'-'))==2,fileParts,'Uni',0)))};
    % Find a subject in this (assumed format: e.g. EB001)
    thisSubj = fileParts{find(cell2mat(cellfun(@(X) length(X)==5,fileParts,'Uni',0)))};

    EphysServerFolder = fullfile(ServerDir,thisSubj,thisDate,'ephys');
    if ~exist('EphysServerFolder')
        mkdir(EphysServerFolder) % Create ephys directory on server if not yet exists
    end
    
    % Check folder for existing .cbin files on server
    cbinfiles = dir(fullfile(EphysServerFolder,'**\*ap*.cbin'));
    if ~isempty(cbinfiles) % This is the already exists loop
        if length(cbinfiles)>1
            disp('What''s happening?!')
            keyboard
        end
        disp('Compressed data already exists on the server...')
        % Check the local version
        LocalFile = dir(fullfile(LocalDir,'**',cbinfiles.name));

        Ok2DeleteLocal = 1;
        % Compare bytes
        if LocalFile.bytes ~= cbinfiles.bytes
            Ok2DeleteLocal = 0;
            disp('What''s happening?!')
            keyboard
        end

        % Check how long cbin file is already on server
        if duration(datetime('now')-datetime(cbinfiles.date))<duration(XCloneDays*24,0,0)
            Ok2DeleteLocal = 0;
            disp(['Copied to server less than ' num2str(XCloneDays) ' days, do not yet delete local copy...'])
        end

        if Ok2DeleteLocal % Delete local copy
            keyboard
            disp(['All seems good, delete local copy...'])
            delete(LocalFile)
        end

    else % This is the compress & copy loop
        % For convenience, first extract the sync from the bin file
        disp('Extracting sync file...')
        [Imecmeta] = ReadMeta2(localEphysFiles.folder,'ap');
        nchan = strsplit(Imecmeta.acqApLfSy,',');
        nChansInFile = str2num(nchan{1})+str2num(nchan{3});
        syncDatImec = extractSyncChannel(localEphysFiles.folder, nChansInFile, nChansInFile); %Last channel is sync

        disp('Compressing data...')
        % Now compress .ap.bin to .ap.cbin
        % Use python integration
        success = pyrunfile("MTSComp_From_Matlab.py","success",datapath = strrep(fullfile(localEphysFiles.folder,localEphysFiles.name),'\','/'),...
            JsonPath =  strrep(fullfile(localEphysFiles.folder,strrep(localEphysFiles.name,'bin','ch')),'\','/'), savepath = strrep(fullfile(localEphysFiles.folder,strrep(localEphysFiles.name,'bin','cbin')),'\','/'))

        if ~success
            keyboard
        else
            disp('All compressed... copy to server now')
            % Copy all data in this folder, except the raw .bin data, to
            % the server
            tmpfiles = dir(localEphysFiles.folder);
            tmpfiles(cellfun(@(X) strfind(X,'*ap.bin')))=[]; %remove the raw bin file from this list
            % Copy the rest
            copyfile(tmpfiles,EphysServerFolder)
        end
    end
end

if TurnOffDesktop
    disp('Great, we''re done. Let''s turn of the computer for now...')
    system('shutdown -s')
end
