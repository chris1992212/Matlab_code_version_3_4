%Cal_ncor
%------------------------------------------
%       Load Data
%------------------------------------------
% load E:\Feng\data\Ian\raw_987987\k-space.mat
% load E:\Feng\data\Ian\raw_987654\k-space.mat
load E:\Feng\data\Ian\raw_OLD_CABLE\k-space.mat
[nPE nFE nCH nSL] = size(data);

mtxCor = zeros(nCH,nCH,nSL);
mtxCov = zeros(nCH,nCH,nSL);
for sl = 1:nSL
    Img = zeros(nPE,nFE,nCH);
    for ch=1:nCH
        Img(:,:,ch)=K2Img(data(:,:,ch,sl));
    end
    [ind] = Get_noise_region(Img,0);
    [cv, C] = ncov(Img,ind);
    mtxCor(:,:,sl) = C;
    mtxCov(:,:,sl) = cv;
end

% save res_987987 mtxCor mtxCov
% save res_987654 mtxCor mtxCov
save res_old_cable mtxCor mtxCov