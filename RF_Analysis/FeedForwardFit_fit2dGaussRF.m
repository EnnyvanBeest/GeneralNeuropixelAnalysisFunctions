function [thisparams,thisgaussian,bestbeta,MSE,coeff,pvalcorr] = fit2dGaussRF(EvokedResponse)

% Create dataset with possible Gaussians
nX = size(EvokedResponse,2);
nY = size(EvokedResponse,1);
nCells = size(EvokedResponse,3);
Centroids = combvec(1:nX,1:nY);
RFSizes = [0.2:0.2:max([nX,nY])/3]; %RF sizes in Full-width half-maximum

ngauss = length(Centroids)*length(RFSizes);
% stimmask = zeros(nX*nY,nX*nY);
% for i = 1:nX*nY
%     stimmask(i,i) = 1;
% end

%% Make all Gaussians - define center with fixed RF Size
disp(['Creating ' num2str(ngauss) ' Gaussians...'])
Gaussians = nan(nY,nX,ngauss);
Gausparams = nan(3,ngauss); %xcenter,Ycenter,FWHM
nrfsizes = length(RFSizes);
countid = 1;
for centerid = 1:size(Centroids,2)
    tmpgauss = nan(nY,nX,nrfsizes);
    tmpparam = nan(3,nrfsizes);
    parfor szid=1:nrfsizes
        gw = normpdf(1:nX,Centroids(1,centerid),RFSizes(szid));
        gv = normpdf(1:nY,Centroids(2,centerid),RFSizes(szid));
        Z = gv'*gw;
        Z = Z./sum(sum(Z)); %Could normalise by sum
        tmpgauss(:,:,szid) = Z;
        tmpparam(:,szid) = [Centroids(:,centerid); RFSizes(szid)];
    end
    Gaussians(:,:,countid:countid+nrfsizes-1) = tmpgauss;
    Gausparams(:,countid:countid+nrfsizes-1) = tmpparam;
    countid = countid+nrfsizes;
end

%% Find Gaussian with largest weight

disp(['Finding best fit for all ' num2str(nCells) ' units...'])

options = optimset('fminsearch');
options.MaxFunEvals = 1000;
beta0 = 1;
for nid = 1:nCells
    parfor gid = 1:ngauss
        tmpfunc = @(X) nanmean((reshape(EvokedResponse(:,:,nid)-abs(X)*Gaussians(:,:,gid),nY*nX,[])).^2,1);
        [beta(gid,nid),mse(gid,nid)]=fminsearch(tmpfunc,beta0,options);
    end
end
%Save parameters
disp('Saving out parameters...')

[MSE,minid] = nanmin(mse,[],1); % Mininmal ME is the one you want
bestbeta = abs(beta(minid));
thisparams = Gausparams(:,minid);
thisgaussian = Gaussians(:,:,minid);

%statistics
ypred = thisgaussian.*permute(repmat(bestbeta',[1,nY,nX]),[2,3,1]);
ypred = reshape(ypred,nY*nX,nCells);
EvokedResponse = reshape(EvokedResponse,nY*nX,nCells);
[coeff,pvalcorr] = arrayfun(@(X) corr(ypred(:,X),EvokedResponse(:,X)),1:nCells,'UniformOutput',0);
coeff = cell2mat(coeff);
pvalcorr = cell2mat(pvalcorr);

end