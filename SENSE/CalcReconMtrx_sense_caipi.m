function [RecMtx,g_map] = CalcReconMtrx_sense_caipi(senseMap,RegI,wrap_m,Ksi,Par,flag_CAIPI)
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
    [rnPE nPE] = size(wrap_m{1});%reduced number of PE
    [nSL rnSL] = size(wrap_m{2});%reduced number of PE
    
    [nPE nSL nFE nCH] = size(senseMap);
    invKsi = inv(Ksi);
    RecMtx = zeros(nPE,nSL,nFE,nCH);
    g_map = ones(nPE,nSL,nFE);
    
    Rf_y = Par.Rf_y;
    Rf_z = Par.Rf_z;
    Rf = Par.Rf;
    for SL = 1:rnSL
        if(flag_CAIPI==0)
            % ************** R=X - SENSE ************** 
            wSL = find(wrap_m{2}(:,SL));
        elseif(flag_CAIPI==1)
            % ************** R=X - CAIPI ************** 
            switch(Rf)
                case 2
                    wSL_1 = find(wrap_m{2}(:,SL));
                    if SL>rnSL/Rf
                        wSL_2 = circshift(wSL_1,[1 0]);
                    else
                        wSL_2 = wSL_1;
                    end
                case 4
                    wSL_1 = repmat(find(wrap_m{2}(:,SL)),[Rf_z 1]);
                    if SL>rnSL/Rf_z
                        wSL_2 = wSL_1;
                    else
                        wSL_2 = circshift(wSL_1,[1 0]);
                    end
            end
        end
        for PE = 1:rnPE
            wPE = find(wrap_m{1}(PE,:));
            if(flag_CAIPI==0)
                % ************** R=X - SENSE **************                
                Rf = max(size(wPE))*max(size(wSL));
            elseif(flag_CAIPI==1)
                % ************** R=X - CAIPI **************
                if PE>rnPE/Rf_y
                    wSL = circshift(wSL_2,[1 0]);
                else
                    wSL = wSL_2;
                end                        
            end             
            for FE = 1:nFE
                if(flag_CAIPI==0)
                    % ************** R=X - SENSE **************  
                    SH = conj(squeeze(senseMap(wPE,wSL,FE,:)));
                    SH = reshape(SH,Rf,nCH);
                    R = diag(reshape(RegI(wPE,wSL,FE),Rf,1));
                elseif(flag_CAIPI==1)
                    % ************** R=X - CAIPI **************  
                    SH = [];
                    for rf=1:Rf
                        SH = [SH;reshape(conj(squeeze(senseMap(wPE(rf),wSL(rf),FE,:))),[1 nCH])];
                    end
                    SH = reshape(SH,Rf,nCH);
                    regI=[];
                    for rf=1:Rf
                        regI=[regI;RegI(wPE(rf),wSL(rf),FE)];
                    end
                    R = diag(reshape(regI,Rf,1));
                end
%                 if Rf ==1
%                     SH = SH.';
%                 end
%                 sprintf('PE=%d, FE=%d',PE,FE)
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
                if(flag_CAIPI==0)
                    % ************** R=X - SENSE **************
                    RecMtx(wPE,wSL,FE,:)= reshape(C,max(size(wPE)),max(size(wSL)),nCH);
                    g_map(wPE,wSL,FE)= reshape(g_tmp,max(size(wPE)),max(size(wSL)));
                elseif(flag_CAIPI==1)
                    % ************** R=X - CAIPI **************
                    for rf=1:Rf
                        RecMtx(wPE(rf),wSL(rf),FE,:) = C(rf,:);
                        g_map(wPE(rf),wSL(rf),FE) = g_tmp(rf);
                    end
                end                
                clear A B C SH R ind g_tmp
            end %for FE
        end %for PE
    end %for SL
end

return