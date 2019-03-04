% This script is used to get the multi-channel Data for reconstruction from original
% acquisition Data
% Input: imgDataCopy_FA1, imgDataCopy_FA1,PCICCMap_mpr, QBC_mpr;
% Output: raw3d_CoilImg,sense_map3d,body_coil3d,CSM;
% Modified by Miao on 14 Feb 2019
function [raw3d_CoilImg, sense_map3d, body_coil3d, CSM] = ...
    Multi_channel_address(imgDataCopy_FA1, imgDataCopy_FA2,PCICCMap_mpr, QBC_mpr)
[nFE, nPE, nmSL, nCH] = size(imgDataCopy_FA1);
% DATA preprocessing for 16 or 24 channels. The original data is uniformly set to 16 channels.
if nCH > 16
   imgDataCopy_FA1 = imgDataCopy_FA1(:,:,:,1:16);
   imgDataCopy_FA2 = imgDataCopy_FA2(:,:,:,1:16);
   PCICCMap_mpr = PCICCMap_mpr(:,:,:,1:16);
   nCH = 16;
end
% split the data
nEcho = 3;
nSL = nmSL/nEcho;
if mod(nSL,nEcho)
    error('The original data is the combination of SL and 3 echo for FA1 and FA2!');
    return;
end
% [nFE, nPE, nSL, nCH] = size(PCICCMap_mpr);
% [nFE, nPE, nSL] = size(QBC_mpr);
raw3d_CoilImg = zeros(nFE, nPE, nSL, nCH, nEcho*2,'single');
for echo = 1:3
    raw3d_CoilImg(:,:,:,:,echo) = single(imgDataCopy_FA1(:,:,1+(echo-1)*48:48*echo,:));
end
clear imgDataCopy_FA1
for echo = 4:6
    raw3d_CoilImg(:,:,:,:,echo) = single(imgDataCopy_FA2(:,:,1+(echo-4)*48:48*(echo-3),:));
end
clear imgDataCopy_FA2
% [nFE_ori,nPE_ori,nSL,nCH,nEcho] = size(raw3d_CoilImg);
PCICCMap_mpr = PCICCMap_mpr(1:76,1:64,:,:);%nFE,nPE,nSL,nCH
QBC_mpr = QBC_mpr(1:76,1:64,:);
raw3d_CoilImg = permute(raw3d_CoilImg,[2 1 3 4 5]); %nPE,nFE,nSL,nCH,nEcho
sense_map3d = permute(PCICCMap_mpr,[2 3 1 4]); %[1 3 2 4],nPE,nSL,nFE,nCH
body_coil3d = permute(QBC_mpr,[2 3 1]);  %[1 3 2] nPE,nSL,nFE
% CSM = permute(PCICCMap_mpr,[2 3 1 4]);%[rnPE, rnSL, rnFE, nCH]
% CSM = Interp_RecMtx(CSM,324,384,48);%[rnPE, rnSL, rnFE, nCH]
% CSM = crop(CSM,[288 48 384 8]);
% CSM = permute(CSM,[1 3 2 4]);
CSM=0;
end
