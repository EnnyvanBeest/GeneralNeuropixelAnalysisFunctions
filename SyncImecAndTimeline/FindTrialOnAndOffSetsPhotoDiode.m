function [starttrialidx,endtrialidx,trialid] = FindTrialOnAndOffSetsPhotoDiode(Actualtime,PhotodiodeSignal,TrialDurations,ITITimes,TrialOnsetTimes)
% Inputs:
% Timeline
% Raw Photodiode signal, convert to binary
% TrialDurations (length 1xntrials)
% ITI times (length 1xntrials, ITI AFTER to stimulus onset of current
% trial)

tmp = detrend(PhotodiodeSignal);
tmp = (tmp-nanmin(tmp(:)))./(nanmax(tmp(:))-nanmin(tmp(:)));
tresh=0.5;
tmp(tmp>=tresh)=1;
tmp(tmp<tresh) = 0;

% Convert to 1 for each flip (so dark and light square change is 1)
tmpdiff = [0 abs(diff(tmp))'];
tresh = 0.05;
tmpdiff(tmpdiff>=tresh)=1;
tmpdiff(tmpdiff~=1)=0;

% If ITI times are unknown, take at least 10 flips
if nargin<4
    ITITimes = repmat(0.4,size(TrialDurations));
end

if nargin<5
    TrialOnsetEstimates = 0;
else
    TrialOnsetEstimates = 1;
end

% sampling rate of time line
tmSR = 1./nanmedian(diff(Actualtime));
WindowSize = round(0.8*(nanmin(ITITimes).*tmSR));

% Upperlimit for refresh rate frames measured
refreshrateInFrames = quantile((diff(find(tmpdiff==1))),0.5);

% Sliding window with average ITI times to find them approximately
filtereddat = filter((1/WindowSize)*ones(1,WindowSize),1,tmpdiff);
diffTime = [nan, diff(filtereddat)];

figure; ax1 = subplot(3,1,1); plot(Actualtime,tmpdiff); title('Photodiode'); hold on;
ax2 = subplot(3,1,2); plot(Actualtime,filtereddat); title('Filtered Photodiode'); hold on
ax3 = subplot(3,1,3); plot(Actualtime,diffTime); title('Change in Filtered Photodiode'); hold on
linkaxes([ax1,ax2,ax3],'x')

starttrialidx = nan(1,length(TrialDurations));
currentlyat = 1;
endtrialidx = nan(1,length(TrialDurations));
TrialDurationsTL = nan(1,length(TrialDurations));
ITIDurationsTL = nan(1,length(TrialDurations));
sessionstart = find(tmpdiff==1,1,'first');
for trialid = 1:length(TrialDurations)
    
    if isempty(find(diffTime(currentlyat+1:end)>0,1,'first')+currentlyat)
        disp('Cannot find anymore trial onsets, perhaps timeline crashed...')
        disp('Continue for now with what we can find...')
        trialid = trialid-1;
        break
    end
    if TrialOnsetEstimates && trialid<length(TrialDurations)
        tmp = find(diffTime(currentlyat+round(0.05*ITITimes(trialid)*tmSR):sessionstart+round(TrialOnsetTimes(trialid+1)*tmSR))>0,1,'first')+round(0.5*ITITimes(trialid)*tmSR)+currentlyat-1;
        if isempty(tmp)
            continue
        else
            starttrialidx(trialid) = tmp;
        end
    else
        if trialid==1 %not always a full ITI recorded
            starttrialidx(trialid) = find(diffTime(currentlyat:end)>0,1,'first')+currentlyat-1;
        else
            starttrialidx(trialid) = find(diffTime(currentlyat+round(ITITimes(trialid)*tmSR):end)>0,1,'first')+currentlyat+round(ITITimes(trialid)*tmSR)-1;
        end
    end
    if trialid==1
        ITIDurationsTL(trialid) = (starttrialidx(trialid)-sessionstart)./tmSR;
    else
        ITIDurationsTL(trialid) = (starttrialidx(trialid)-currentlyat)./tmSR;
    end
        
    
    currentlyat =  starttrialidx(trialid);
    if TrialOnsetEstimates && trialid<length(TrialDurations)
        tmp = find(filtereddat(currentlyat+round(0.1*TrialDurations(trialid)*tmSR):sessionstart+round(TrialOnsetTimes(trialid+1)*tmSR))==0,1,'first')+currentlyat+round(0.1*TrialDurations(trialid)*tmSR)-1-WindowSize;
        if isempty(tmp)
            %Try different method
            tmpsig = tmpdiff(currentlyat:sessionstart+round(TrialOnsetTimes(trialid+1)*tmSR));
            SmallerWindowSize = round(0.2.*tmSR);
            filteredtmp = filter((1/SmallerWindowSize)*ones(1,SmallerWindowSize),1,tmpdiff);            
            tmp = find(filteredtmp(currentlyat+round(0.1*TrialDurations(trialid)*tmSR):sessionstart+round(TrialOnsetTimes(trialid+1)*tmSR))==0,1,'first')+currentlyat+round(0.1*TrialDurations(trialid)*tmSR)-1-SmallerWindowSize;
            if isempty(tmp)
                continue
            else
                endtrialidx(trialid)= tmp;                
            end            
        else
            endtrialidx(trialid)= tmp;
        end
    else
        endtrialidx(trialid) = find(filtereddat(currentlyat:end)==0,1,'first')+currentlyat-1-WindowSize;
    end
    currentlyat =  endtrialidx(trialid);
    
    TrialDurationsTL(trialid) = (endtrialidx(trialid)-starttrialidx(trialid))./tmSR;
    for subid=1:3
        subplot(3,1,subid)
        if currentlyat-50000>0 && currentlyat+25000<length(Actualtime)
            xlim([Actualtime(currentlyat-50000) Actualtime(currentlyat+25000)])
        elseif currentlyat-50000>0
            xlim([Actualtime(currentlyat-50000) Actualtime(end)])
        else
            try
                xlim([Actualtime(1) Actualtime(currentlyat+25000)])
            catch
                xlim([Actualtime(1) Actualtime(end)])
            end
        end
        line([Actualtime(starttrialidx(trialid)) Actualtime(starttrialidx(trialid))],get(gca,'ylim'),'color',[0 1 0],'LineWidth',2)
        line([Actualtime(endtrialidx(trialid))  Actualtime(endtrialidx(trialid))],get(gca,'ylim'),'color',[1 0 0],'LineWidth',2)
    end
    drawnow
end

for subid=1:3
    subplot(3,1,subid)
    xlim([Actualtime(1) Actualtime(end)])
end

drawnow
ITIDurationsTL(1) = 0;
if trialid<length(TrialDurations)
    %% Do trialdurations match?
    [Cdur,Lagdur] = xcorr(TrialDurations',TrialDurationsTL(1:trialid));
    [CIti,LagIti] = xcorr(ITITimes,ITIDurationsTL(1:trialid));
    
    
    [maxval,maxid1] = max(Cdur); Lagdur(maxid1)
    [maxval,maxid2] = max(CIti); Lagdur(maxid2)
    [maxval,maxid3] = max(CIti+Cdur*2); Lagdur(maxid3)
    startedontrial = mode([Lagdur(maxid1) Lagdur(maxid2) Lagdur(maxid3)]);
    
    figure;
    subplot(3,1,1); plot(Lagdur,Cdur); title('CrossCorrelation Trial Durations'); hold on; line([startedontrial, startedontrial],get(gca,'ylim'),'color',[0 0 0])
    subplot(3,1,2); plot(LagIti,CIti); title('CrossCorrelation  ITIs'); hold on; line([startedontrial, startedontrial],get(gca,'ylim'),'color',[0 0 0])
    subplot(3,1,3); plot(LagIti,CIti+Cdur); title('CrossCorrelation  Sum'); hold on; line([startedontrial, startedontrial],get(gca,'ylim'),'color',[0 0 0])
    
    starttrialidxnew = nan(1,length(TrialDurations));
    starttrialidxnew(startedontrial+1:startedontrial+trialid) = starttrialidx(1:trialid);
    endtrialidxnew = nan(1,length(TrialDurations));
    endtrialidxnew(startedontrial+1:startedontrial+trialid) = endtrialidx(1:trialid);
    
    starttrialidx=starttrialidxnew;
    endtrialidx = endtrialidxnew;
    
    AcceptableYN = '';
    while ~ismember(AcceptableYN,{'y','n'})
        drawnow
        AcceptableYN = input('Is this alignment acceptable? (y/n)','s');
        if strcmpi(AcceptableYN,'n')
            disp('Not acceptable, skipping session')
            starttrialidx=[];
            endtrialidx = [];
            break
        elseif strcmp(AcceptableYN,'y')
            disp('Acceptable, continuing...')
        end
    end
    
end