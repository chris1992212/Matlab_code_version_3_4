function [RecMtx,g_map] = CalcReconMtrx(senseMap,RegI,wrap_m,Ksi)
%This function is used to calculation reconstruction matrix using
%regularization matrix 
%INPUT
%   senseMap: sensitivity maps [nPE nFE nCH]
%   RegI: regularization image [nPE nFE], 1./R
%   wrap_m: SENSE warping map [rnPE nPE]
%   Ksi: noise correlation matrix [nCH nCH]
%OUTPUT
%   RecMtx: reconstruction matrix [nPE nFE nCH]
%Reference
%   RS: Reconstruction Page 235~236
[rnPE nPE] = size(wrap_m);%reduced number of PE
[nPE nFE nCH] = size(senseMap);
invKsi = inv(Ksi);
RecMtx = zeros(nPE,nFE,nCH);
g_map = ones(nPE,nFE);
for PE = 1:rnPE
    wPE = find(wrap_m(PE,:));
    for FE = 1:nFE
        SH = conj(squeeze(senseMap(wPE,FE,:)));
        if size(wPE)==[1 1]
            SH = SH.';
        end
        R = diag(RegI(wPE,FE));
        [SH,R,C,ind]=RmvZero(SH, R); 
        g_tmp = ones(size(C,1),1);
        if ~isempty(SH)
            B = SH*invKsi;
            A = B*SH' + R;
            inv_A = inv(A);
            C(ind,:) = inv_A*B;
            g_tmp(ind)=sqrt(diag(inv_A).*diag(A));
        end
        RecMtx(wPE,FE,:)=C;
        g_map(wPE,FE)=g_tmp;
        clear A B C SH R ind g_tmp
    end %for FE
end %for PE
% [RecMtx] = smoothMtx(RecMtx);


return