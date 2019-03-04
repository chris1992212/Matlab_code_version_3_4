%demo_CS_SENSE

%-------------------------------------------------------------------
% Define parameters
%-------------------------------------------------------------------
[Par]=define_Par();


%-------------------------------------------------------------------
% load Data
%-------------------------------------------------------------------
nsl = 5;
load G:\FengBackup\Data\philips\Gainesville\3T\SFSENSE\raw_029\ImgKsl1_8.mat
ImgK = squeeze(ImgK(:,:,:,nsl));
%segment FE 
Kdata = flipdim(permute(ImgK,[2 1 3]),1);
[nPE, nFE, nCH] = size(Kdata);
clear ImgK
load G:\FengBackup\Data\philips\Gainesville\3T\SFSENSE\raw_029\Phi.mat
%-------------------------------------------------------------------
%---Caculate sensitivity maps---
%-------------------------------------------------------------------
%load low resolution bodycoil
load G:\FengBackup\Data\philips\Gainesville\3T\SFSENSE\raw_027\ImgK.mat
bodyK = squeeze(data(:,:,nsl));
clear data
[bodyK]=flipdim(permute(bodyK,[2 1]),1);
bodycoil = K2Img(bodyK);
clear bodyK
%load low resolution multi-channel data
load G:\FengBackup\Data\philips\Gainesville\3T\SFSENSE\raw_026\ImgK.mat
LresK = squeeze(data(:,:,:,nsl));
clear data
[LresK]=flipdim(permute(LresK,[2 1 3]),1);
for ch = 1:nCH
    LresImg(:,:,ch) = K2Img(LresK(:,:,ch));
end
%calculate sensitivity maps
LSS = sqrt(sum(abs(LresImg).^2,3));
[ind_big,Im]=sep_1d(abs(LSS(:)),1);
mask = zeros(size(LSS));
mask(ind_big{1})=1;
% filled_mask = imfill(mask,'holes');
se = strel('disk',2); %default is 2
filled_mask = imfill(medfilt2(imdilate(double(mask),se)),'holes');
SSM = CalcCoilSensitivities(LresImg,bodycoil,filled_mask);

% [SSM ] = Interp_RecMtx(SSM ,nPE,nFE);
% [bodycoil] = Interp_RecMtx(bodycoil,nPE,nFE);

clear LresK

%-------------------------------------------------------------------
% Initial Reconstruction
%-------------------------------------------------------------------
ParK = zeros(nPE,nFE,nCH);
ParK(1:Par.Rf:nPE,:,:) = Kdata(1:Par.Rf:nPE,:,:);

tic
[CompI,RecMtx,g_map]=RegPreSENSE(ParK,Phi,SSM, bodycoil,Par);
toc
% CompI = CompI.*Par.Rf;
%-------------------------------------------------------------------
% self-feeding Sparse SENSE
%-------------------------------------------------------------------
% tic;[Inew,senseMap]=Self_Sparse_SENSE(CompI,ParK,SSM,g_map,Phi,Par);toc
% addpath('E:\projects\ktGROWL\files\GROWL');
tic;[Inew,senseMap]=Sparse_GoSENSE(CompI,ParK,SSM,g_map,Phi,Par);toc
%-------------------------------------------------------------------
% Show results
%-------------------------------------------------------------------
real_I = zeros(nPE,nFE);
for ch=1:nCH
    real_I = real_I + K2Img(Kdata(:,:,ch)).*conj(senseMap(:,:,ch));
end
figure;imagesc(abs(Inew)',[0 3]);axis off;axis equal;colormap(gray)
% figure;imagesc(abs(CompI)',[0 3]);axis off;axis equal;colormap(gray)
% figure;imagesc(abs(abs(Inew)-abs(real_I))',[0 0.6]);axis off;axis equal;colormap(gray)
% figure;imagesc(abs(T(161:224,201:264)),[0 3]);axis off;axis equal;colormap(gray)
norm(abs(Inew)-abs(real_I),'fro')/norm(real_I,'fro')
norm(abs(CompI)-abs(real_I),'fro')/norm(real_I,'fro')
