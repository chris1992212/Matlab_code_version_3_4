function [Inew,senseMap]=Self_Sparse_SENSE(CompI,ParK,SSM,g_map,Phi,Par)
%INPUTS
%   CompI: Initial reconstruction by on-line SENSE, [nPE nFE]
%   ParK: Partially acquired k-space [nPE nFE nCH]
%   senseMap: sensitivity maps from system [nsPE nsFE nCH]
%   g_map: g_factor map from system [nPE nFE ]
%   Phi: noise correlation matrix [nCH nCH]
%   Par: parameters 
%OUTPUTS
%   Inew: Reconstructed image
%-------------------------------------------------------------------
%  Feng Huang, Aug 12, 2010
%-------------------------------------------------------------------
if ndims(CompI)==2
    [nPE nFE ] = size(CompI);
    %---find the background-----
    [ind_big,Im]=sep_1d(abs(CompI(:)),1);
    mask = zeros(nPE,nFE);
    mask(ind_big{1})=1;
    [senseMap] = Interp_RecMtx(SSM ,nPE,nFE);
elseif ndims(CompI)==3
    [nPE nSL nFE ] = size(CompI);
    %---find the background-----
    [ind_big,Im]=sep_1d(abs(CompI(:)),1);
    mask = zeros(nPE,nSL,nFE);
    mask(ind_big{1})=1;
    [senseMap] = Interp_RecMtx(SSM ,nPE,nFE,nSL);
%     SS = sqrt(sum(abs(senseMap).^2,4));
%     for ch = 1:size(senseMap,4)
%         tmpSM = zeros(nPE,nFE,nSL);
%         tmpImg = senseMap(:,:,:,ch);
%         tmpSM = tmpImg./SS;
%         senseMap(:,:,:,ch) = tmpSM;
%     end
%     senseMap(isnan(senseMap))=0;
end
if (mean(abs(g_map(find(mask))))>Par.Lg)
    
    %Smask = zeros(nPE,nFE);
    %Smask(find(abs(g_map)>Par.Lg)) = 1;
%     Par.aTV = Par.aTV*(mean(abs(g_map(find(mask))))); %parameters for TV term
    Par.weight = abs(abs(g_map)-1)*(mean(abs(g_map(find(mask)))));
    %-------------------------------------------------------------------
    % Image denoising
    %-------------------------------------------------------------------
    disp('denoising')
    tic
    if ndims(CompI)==2
        [nI,Out] = TVDenoise(CompI,Par);
    elseif ndims(CompI)==3
        [nI,Out] = TVDenoise3D(CompI,Par);
%         nI = zeros(size(CompI));
%         sPar = Par;
%         for FE = 1:nFE
%             sPar.weight = squeeze(Par.weight(:,:,FE));
%             CompI2D = squeeze(CompI(:,:,FE));
%             if max(abs(CompI2D))==0
%                 nI(:,:,FE) = CompI2D;
%             else
%                 [nI(:,:,FE),Out] = TVDenoise(CompI2D,sPar);
%             end
%         end
    end 
    toc
    disp('denoising')
    % %-------------------------------------------------------------------
    % Prior information regularized SENSE
    % %------------------------------------------------------------------- 
    if norm(abs(nI(:))-abs(CompI(:)),'fro')/norm(abs(CompI(:)),'fro') > 0% 0.02
         tic;[Inew,RecMtx]=PriRegSENSEtmp(ParK,nI,senseMap,SSM,Phi,Par);toc
%         [Inew]=PriRegSENSE(ParK,nI,senseMap,Phi,Par);
%         if ~isfield(Par,'innitr')
%             Par.innitr = 15; %number of iternal iterations
%         end
%         if ~isfield(Par,'eTol')
%             Par.eTol = 1e-4; %determine if it is converged
%         end
%         [Inew]=PriRegSENSE_LSBB(ParK,nI,senseMap,Phi,Par);
    else
        Inew = CompI;
    end
    %Inew = Inew.*Smask + CompI.*(1-Smask);
    %     Inew = nI;
else
    Inew = CompI;
end
