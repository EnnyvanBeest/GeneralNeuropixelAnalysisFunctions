function [starttrialidx,endtrialidx,trialid] = FindTrialOnAndOffSetsPhotoDiode(Actualtime,PhotodiodeSignal,TrialDurations,ITITimes,TrialOnsetTimes)
showplot = 0;
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

tmpdif = [0 diff(filtereddat==0)];
starttrialidx = find(tmpdif==-1);
endtrialidx = cell2mat(arrayfun(@(X) find(filtereddat(X:end)==0,1,'first')+X-1-WindowSize,starttrialidx,'UniformOutput',0));
if length(endtrialidx)<length(starttrialidx)
    endtrialidx(length(endtrialidx)+1:length(starttrialidx)) = length(Actualtime);
end
TrialDurationsTL = Actualtime(endtrialidx)-Actualtime(starttrialidx);

for trialid = 1:length(starttrialidx)
    if trialid==1
        ITIDurationsTL(trialid) = (starttrialidx(trialid))./tmSR;
    else
        ITIDurationsTL(trialid) = (starttrialidx(trialid)-endtrialidx(trialid-1))./tmSR;
    end
    for subid=1:3
        subplot(3,1,subid)
        if starttrialidx(trialid)-50000>0 && starttrialidx(trialid)+25000<length(Actualtime)
            xlim([Actualtime(starttrialidx(trialid)-50000) Actualtime(starttrialidx(trialid)+25000)])
        elseif starttrialidx(trialid)-50000>0
            xlim([Actualtime(starttrialidx(trialid)-50000) Actualtime(end)])
        else
            try
                xlim([Actualtime(1) Actualtime(starttrialidx(trialid)+25000)])
            catch
                xlim([Actualtime(1) Actualtime(end)])
            end
        end
        line([Actualtime(starttrialidx(trialid)) Actualtime(starttrialidx(trialid))],get(gca,'ylim'),'color',[0 1 0],'LineWidth',2)
        line([Actualtime(endtrialidx(trialid))  Actualtime(endtrialidx(trialid))],get(gca,'ylim'),'color',[1 0 0],'LineWidth',2)
    end
    if showplot
        
        drawnow
    end
end

%
% starttrialidx = nan(1,length(TrialDurations));
% currentlyat = 1;
% endtrialidx = nan(1,length(TrialDurations));
% TrialDurationsTL = nan(1,length(TrialDurations));
% ITIDurationsTL = nan(1,length(TrialDurations));
% sessionstart = find(tmpdiff==1,1,'first');
% for trialid = 1:length(TrialDurations)
%     if isempty(find(diffTime(currentlyat+1:end)>0,1,'first')+currentlyat)
%         disp('Cannot find anymore trial onsets, perhaps timeline crashed...')
%         disp('Continue for now with what we can find...')
%         trialid = trialid-1;
%         break
%     end
%     if TrialOnsetEstimates && trialid<length(TrialDurations)
%         tmp = find(diffTime(currentlyat+round(0.05*ITITimes(trialid)*tmSR):sessionstart+round(TrialOnsetTimes(trialid+1)*tmSR))>0,1,'first')+round(0.5*ITITimes(trialid)*tmSR)+currentlyat-1;
%         if isempty(tmp)
%             continue
%         else
%             starttrialidx(trialid) = tmp;
%         end
%     elseif trialid==1
%         starttrialidx(trialid) = find(diffTime(currentlyat:end)>0,1,'first')+currentlyat-1;
%         if ~any(diffTime(starttrialidx(trialid)+1:starttrialidx(trialid)+tmSR.*3/4)>0) % Sometimes one flip happens before actual trial onset --> skip this frame (assuming there's at least some >1 refresh rate
%             starttrialidx(trialid) = find(diffTime(starttrialidx(trialid)+1:end)>0,1,'first')+starttrialidx(trialid);
%         end
%         TrialOnsetTimes(trialid)=Actualtime(starttrialidx(trialid));
%     else
%         starttrialidx(trialid) = find(diffTime(currentlyat+round(ITITimes(trialid)*tmSR):end)>0,1,'first')+currentlyat+round(ITITimes(trialid)*tmSR)-1;
%     end
%     if trialid==1
%         ITIDurationsTL(trialid) = (starttrialidx(trialid)-sessionstart)./tmSR;
%     else
%         ITIDurationsTL(trialid) = (starttrialidx(trialid)-currentlyat)./tmSR;
%     end
%
%
%     currentlyat =  starttrialidx(trialid);
%     if TrialOnsetEstimates && trialid<length(TrialDurations)
%         tmp = find(filtereddat(currentlyat+round(0.1*TrialDurations(trialid)*tmSR):sessionstart+round(TrialOnsetTimes(trialid+1)*tmSR))==0,1,'first')+currentlyat+round(0.1*TrialDurations(trialid)*tmSR)-1-WindowSize;
%         if isempty(tmp)
%             %Try different method
%             tmpsig = tmpdiff(currentlyat:sessionstart+round(TrialOnsetTimes(trialid+1)*tmSR));
%             SmallerWindowSize = round(0.2.*tmSR);
%             filteredtmp = filter((1/SmallerWindowSize)*ones(1,SmallerWindowSize),1,tmpdiff);
%             tmp = find(filteredtmp(currentlyat+round(0.1*TrialDurations(trialid)*tmSR):sessionstart+round(TrialOnsetTimes(trialid+1)*tmSR))==0,1,'first')+currentlyat+round(0.1*TrialDurations(trialid)*tmSR)-1-SmallerWindowSize;
%             if isempty(tmp)
%                 continue
%             else
%                 endtrialidx(trialid)= tmp;
%             end
%         else
%             endtrialidx(trialid)= tmp;
%         end
%     else
%         endtrialidx(trialid) = find(filtereddat(currentlyat:end)==0,1,'first')+currentlyat-1-WindowSize;
%     end
%     currentlyat =  endtrialidx(trialid);
%
%     TrialDurationsTL(trialid) = (endtrialidx(trialid)-starttrialidx(trialid))./tmSR;
%
%     for subid=1:3
%         subplot(3,1,subid)
%         if currentlyat-50000>0 && currentlyat+25000<length(Actualtime)
%             xlim([Actualtime(currentlyat-50000) Actualtime(currentlyat+25000)])
%         elseif currentlyat-50000>0
%             xlim([Actualtime(currentlyat-50000) Actualtime(end)])
%         else
%             try
%                 xlim([Actualtime(1) Actualtime(currentlyat+25000)])
%             catch
%                 xlim([Actualtime(1) Actualtime(end)])
%             end
%         end
%         line([Actualtime(starttrialidx(trialid)) Actualtime(starttrialidx(trialid))],get(gca,'ylim'),'color',[0 1 0],'LineWidth',2)
%         line([Actualtime(endtrialidx(trialid))  Actualtime(endtrialidx(trialid))],get(gca,'ylim'),'color',[1 0 0],'LineWidth',2)
%     end
%     drawnow
% end

for subid=1:3
    subplot(3,1,subid)
    xlim([Actualtime(1) Actualtime(end)])
end
drawnow

%% Try to align Trial Durations
ntrials = min([length(TrialDurations) length(TrialDurationsTL)]);
if length(TrialDurations)~=length(TrialDurationsTL)
    display(['Cannot find the appropriate number of trials...'])
    
    [Cdur,Lagdur] = xcorr(TrialDurations',TrialDurationsTL);
    [CIti,LagIti] = xcorr(ITITimes,ITIDurationsTL);
    
    [maxval,maxid1] = max(Cdur); Lagdur(maxid1)
    [maxval,maxid2] = max(CIti); Lagdur(maxid2)
    [maxval,maxid3] = max(CIti+Cdur*2); Lagdur(maxid3)
    startedontrial = mode([Lagdur(maxid1) Lagdur(maxid2) Lagdur(maxid3)]);
    
    figure;
    subplot(3,1,1); plot(Lagdur,Cdur); title('CrossCorrelation Trial Durations'); hold on; line([startedontrial, startedontrial],get(gca,'ylim'),'color',[0 0 0])
    subplot(3,1,2); plot(LagIti,CIti); title('CrossCorrelation  ITIs'); hold on; line([startedontrial, startedontrial],get(gca,'ylim'),'color',[0 0 0])
    subplot(3,1,3); plot(LagIti,CIti+Cdur); title('CrossCorrelation  Sum'); hold on; line([startedontrial, startedontrial],get(gca,'ylim'),'color',[0 0 0])
    
    starttrialidxnew = nan(1,length(TrialDurations));
    endtrialidxnew = nan(1,length(TrialDurations));
    
    if startedontrial>=0
        starttrialidxnew(startedontrial+1:startedontrial+ntrials) = starttrialidx(1:ntrials);
        endtrialidxnew(startedontrial+1:startedontrial+ntrials) = endtrialidx(1:ntrials);
        TrialDurationsTL = TrialDurationsTL(1:ntrials);
        
    else
        starttrialidxnew(1:ntrials) = starttrialidx(-startedontrial+1:ntrials-startedontrial);
        endtrialidxnew(1:ntrials) = endtrialidx(-startedontrial+1:ntrials-startedontrial);
        TrialDurationsTL = TrialDurationsTL(-startedontrial+1:ntrials-startedontrial);
    end
    starttrialidx=starttrialidxnew;
    endtrialidx = endtrialidxnew;
    
end


%     AcceptableYN = '';
%     while ~ismember(AcceptableYN,{'y','n'})
%         drawnow
%         AcceptableYN = input('Is this alignment acceptable? (y/n)','s');
%         if strcmpi(AcceptableYN,'n')
%             disp('Not acceptable, skipping session')
%             starttrialidx=[];
%             endtrialidx = [];
%             break
%         elseif strcmp(AcceptableYN,'y')
%             disp('Acceptable, continuing...')
%         end
%     end



figure;
scatter(TrialDurations(~isnan(starttrialidx)),TrialDurationsTL)
hold on
line([min(TrialDurations(~isnan(starttrialidx))) max(TrialDurations(~isnan(starttrialidx)))],[min(TrialDurations(~isnan(starttrialidx))) max(TrialDurations(~isnan(starttrialidx)))],'color',[1 0 0])
end
