%demo_CS_SENSE
%-------------------------------------------------------------------
% Locate necessary files
%-------------------------------------------------------------------
addpath('E:\Work\BestProjects\SFSENSE\3DSparseSENSE\files\ImageDenoising\')
addpath('E:\Work\BestProjects\SFSENSE\3DSparseSENSE\files\SENSE\')
addpath('E:\Work\BestProjects\SFSENSE\3DSparseSENSE\files\SensitivityMap\')
addpath('E:\Work\BestProjects\SFSENSE\3DSparseSENSE\files\SparseSENSE\')
%-------------------------------------------------------------------
% Define parameters
%-------------------------------------------------------------------
Rf = 4;
Par.nPE = 240;
Par.nSL = 64;
%parameters for Philips SENSE
% Par.Rf = Rf;
Par.Areg = 1e-6;%0.0001;
Par.alpha = 2;
Par.beta = 1;
Par.Kar = 0;
%parameters for sensitivity maps
Par.flagSM = 0;
%parameters for image denoising
Par.aTV = 0.02;
Par.TVtype = 2;
Par.maxItr = 100;       
Par.gamma = 1.0;        
Par.beta = 10;          
Par.relchg_tol = 5e-4;  
Par.normalize = 1;
%parameters for prior information regularized SENSE
Par.priLmd = 0.1; %weight for regularization
Par.Lg = 1.0; %no change at locations where g-factor is smaller than this value
%-------------------------------------------------------------------
% load Data
%-------------------------------------------------------------------
load C:\Feng\data\3DSENSE\raw_012\ImgK.mat
ParK = permute(data,[2 4 1 3]);
clear data
ParK = ParK(1:4:end,1:2:end,:,:);
[nrPE, nrSL,nFE, nCH] = size(ParK);
load C:\Feng\data\3DSENSE\Phi.mat
%-------------------------------------------------------------------
%---Caculate sensitivity maps---
%-------------------------------------------------------------------
%load low resolution bodycoil
load C:\Feng\data\3DSENSE\raw_014\ImgK.mat
bodyK = permute(data,[2 3 1]);
clear data

bodycoil = K2Img(bodyK);
clear bodyK
%load low resolution multi-channel data
load C:\Feng\data\3DSENSE\raw_013\ImgK.mat
LresK = permute(data,[2 4 1 3]);
clear data

for ch = 1:nCH
    LresImg(:,:,:,ch) = K2Img(LresK(:,:,:,ch));
end
clear LresK
%calculate sensitivity maps
tic
[SSM,bodycoil,filled_mask]=SM_LresImg(LresImg,Par,bodycoil);
toc



%-------------------------------------------------------------------
% Initial Reconstruction
%-------------------------------------------------------------------


tic
[CompI,RecMtx,g_map]=RegPreSENSE(ParK,Phi,SSM, bodycoil,Par);
toc
% figure;imagesc(squeeze(abs(CompI(:,32,:))));colormap(gray)
% figure;imagesc(squeeze(abs(CompI(:,:,120))));colormap(gray)
%-------------------------------------------------------------------
% self-feeding Sparse SENSE
%-------------------------------------------------------------------
[Inew,senseMap]=Self_Sparse_SENSE(CompI,ParK,SSM,g_map,Phi,Par);
% tic;[Inew,senseMap]=Self_Sparse_SENSE(CompI,ParK,SSM,g_map,Phi,Par);toc
% addpath('E:\projects\ktGROWL\files\GROWL');
% tic;[Inew,senseMap]=Sparse_GoSENSE(CompI,ParK,SSM,g_map,Phi,Par);toc
%-------------------------------------------------------------------
% Show results
%-------------------------------------------------------------------

load realI.mat
Inew = Inew*norm(realI(:),'fro')/norm(Inew(:),'fro');
diff = realI - Inew;
err = norm(diff(:),'fro')/norm(realI(:),'fro')
figure;imagesc(squeeze(abs(Inew(61:180,32,61:180)))');colormap(gray);axis off;axis equal
figure;imagesc(squeeze(abs(realI(61:180,32,61:180)))');colormap(gray);axis off;axis equal
figure;imagesc(squeeze(abs(CompI(61:180,32,61:180)))');colormap(gray);axis off;axis equal