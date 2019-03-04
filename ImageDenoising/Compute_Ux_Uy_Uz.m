function [Ux,Uy,Uz]=Compute_Ux_Uy_Uz(U)
%This function is used to calculate finite difference transform along 3 dimenstions
%INPUT
%   U: a 3D image, [nPE nFE nSL]
%OUTPUTs
%   Ux, Uy, Uz: finite difference transform along 3 dimensions

[nPE, nFE, nSL] = size(U);
Ux = zeros(nPE,nFE,nSL);
Ux(:,1:nFE-1,:) = diff(U,1,2);
Ux(:,nFE,:) = U(:,1,:) - U(:,nFE,:);


Uy = zeros(nPE,nFE,nSL);
Uy(1:nPE-1,:,:) = diff(U,1,1);
Uy(nPE,:,:) = U(1,:,:) - U(nPE,:,:);

Uz = zeros(nPE,nFE,nSL);
Uz(:,:,1:nSL-1) = diff(U,1,3);
Uz(:,:,nSL) = U(:,:,1) - U(:,:,nSL);
