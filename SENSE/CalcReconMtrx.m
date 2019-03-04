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
if ~iscell(wrap_m) %1D SENSE
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
else %2D SENSE
    [rnPE , ~] = size(wrap_m{1});%reduced number of PE
    [~, rnSL] = size(wrap_m{2});%reduced number of PE
    
    [nPE, nSL, nFE, nCH] = size(senseMap);
    invKsi = inv(Ksi);
    RecMtx = zeros(nPE,nSL,nFE,nCH,'single');
    g_map = ones(nPE,nSL,nFE,'single');
    
    for SL = 1:rnSL
        wSL = find(wrap_m{2}(:,SL));
        for PE = 1:rnPE
            wPE = find(wrap_m{1}(PE,:));
            Rf = max(size(wPE))*max(size(wSL));
            for FE = 1:nFE
                SH = conj(squeeze(senseMap(wPE,wSL,FE,:)));
                SH = reshape(SH,Rf,nCH);
%                 if Rf ==1
%                     SH = SH.';
%                 end
%                 sprintf('PE=%d, FE=%d',PE,FE)
                R = diag(reshape(RegI(wPE,wSL,FE),Rf,1));
                if size(R,1)==size(SH,1)
                [SH,R,C,ind]=RmvZero(SH, R);
                else
                    pause
                end
                g_tmp = ones(size(C,1),1);
                if ~isempty(SH)
                    B = SH*invKsi;
                    A = B*SH' + R;
                    inv_A = inv(A);
                    C(ind,:) = inv_A*B;
                    g_tmp(ind)=sqrt(diag(inv_A).*diag(A));
                end
                RecMtx(wPE,wSL,FE,:)= reshape(C,max(size(wPE)),max(size(wSL)),nCH);
                g_map(wPE,wSL,FE)= reshape(g_tmp,max(size(wPE)),max(size(wSL)));
%                 clear A B C SH R ind g_tmp
            end %for FE
        end %for PE
    end %for SL
end

return