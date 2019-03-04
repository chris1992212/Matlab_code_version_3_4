function [Inew,RecMtx]=PriRegSENSEtmp_CAIPI(Kdata,PriI,Isens,SSM,Ksi,Par,flag_CAIPI)
%This function is for the prior information regularized SENSE
%This version is only for 3D images
%Image support reduction idea is implemented
%INPUT
%   Kdata: partial k-space data [nacqPE nacqSL nFE nCH]
%   PriI: Prior information of the image, a 2D matrix, [nPE nSL nFE]
%   Isens: sensitivity maps [nPE nSL nFE nCH]
%   Ksi : noise correlation matrix
%   Par: parameters
%       .priLmd: weight Lamdba for the prior informaiton
%OUTPUT
%   Inew: reconstructed image
%-----------------------------------------------------------------
%               Feng Huang, Jan 05, 2011
%      Explicit image support reduction technique is used
%  Aiqi Sun,   Apr 8,  2018 => for CAIPI-SFSENSE recon
%-----------------------------------------------------------------

%-----------------------------------------------------------------
%Produce wrapped image with reduced image support reduction
%-----------------------------------------------------------------
[nacqPE,nacqSL,nFE,nCH]=size(Kdata);

nPE = Par.nPE;
nSL = Par.nSL;

%------ Produce the full FOV wrapped image------ 
while(0)
    
%produce aliasing matrix
[alias_m]=sense_alias_matrix([Par.nPE Par.nSL],[nacqPE nacqSL]);
%produce the full FOV wrapped image
[WrapI] = genWrapI(Kdata,alias_m);

ImgC = zeros(size(Isens));
for ch=1:nCH
    ImgC(:,:,:,ch) = PriI.*Isens(:,:,:,ch);
end
[nWrapI] = genWrapI(ImgC,alias_m);

end

% ************** R=X - CAIPI ************** 
for ch=1:nCH
    WrapI(:,:,:,ch) = K2Img_3d(Kdata(:,:,:,ch)*Par.Rf);
end
ImgC = zeros(size(Isens));
for ch=1:nCH
    ImgC(:,:,:,ch) = PriI.*Isens(:,:,:,ch);
end
Rf_y = 2;
Rf_z = 2;
ParK_00 = fft3c(ImgC);
tmp_KK = zeros(size(ParK_00));
tmp_KK(1:Rf_y*Rf_z:end,1:Rf_z:end,:,:)=ParK_00(1:Rf_y*Rf_z:end,1:Rf_z:end,:,:);
tmp_KK(Rf_y*Rf_z-1:Rf_y*Rf_z:end,2:Rf_z:end,:,:)=ParK_00(Rf_y*Rf_z-1:Rf_y*Rf_z:end,2:Rf_z:end,:,:);
clear ImgC ParK_00;
nWrapI = ifft3c(tmp_KK)*Par.Rf;
clear tmp_KK;
% as(sos(nWrapI,4))
% norm(sos(squeeze(WrapI(:,14,:,:)),3),'fro')
% norm(sos(squeeze(nWrapI(:,14,:,:)*Par.Rf),3),'fro')

for ch=1:nCH
    tmp_nWrapI = nWrapI(:,:,:,ch);
    tmp_WrapI = WrapI(:,:,:,ch);
    nWrapI(:,:,:,ch) = nWrapI(:,:,:,ch)/norm(tmp_nWrapI(:),'fro')*norm(tmp_WrapI(:),'fro');
    nWrapI(:,:,:,ch) = WrapI(:,:,:,ch) - nWrapI(:,:,:,ch);
end

RecMtx = [];%Par.RecMtx;

%-----------------------------------------------------------------
%Produce new reconstruction matrix
%-----------------------------------------------------------------
if isempty(RecMtx)
    %Produce bodycoil for regularization
%     [dx,dy] = gradient(PriI);
    bodycoil = ones(nPE,nSL,nFE);%./(1e-2+Par.weight); %sqrt(dx.^2+dy.^2)
%     bodycoil = Par.weight;%./sqrt(abs(dx).^2+abs(dy).^2);
    %acquired number of PE for lower resolution
    [nrPE,nrSL,nrFE,nCH]=size(SSM);
        
    %acquired number of PE for lower resolution
    ParNrPE = round(nrPE/Par.Rf_y); % round(nacqPE/Par.nPE*nrPE);
    ParNrSL = round(nrSL/Par.Rf_z); % round(nacqSL/Par.nSL*nrSL);

    [wrap_m]=sense_alias_matrix([nrPE nrSL],[ParNrPE ParNrSL]);
    if ParNrPE==nrPE
        wrap_m{1}=eye(nrPE);
    end
    if ParNrSL == nrSL
        wrap_m{2}=eye(nrSL);
    end
    
    if(flag_CAIPI==1)
        switch(Par.Rf)
            case 2
                tmp_wrap = circshift(wrap_m{2},[0 nrSL/Par.Rf]);
                %     figure,imshow(tmp_wrap)
                wrap_m2_orig = wrap_m{2};
                wrap_m{2} = tmp_wrap + eye(nrSL);
                %     figure,imshow(wrap_m{2})
            case 4
                tmp_wrap = repmat(eye(nrPE/Par.Rf_y),[1 2]);
                %     figure,imshow(tmp_wrap)
                wrap_m1_orig = wrap_m{1};
                wrap_m{1} = wrap_m{1} + tmp_wrap;
                %     figure,imshow(wrap_m{1})
        end
    end
    
    [bodycoil] = abs(Interp_RecMtx(bodycoil,nrPE,nrFE,nrSL));
    
    Par.Areg = 100;% 1e-3;%5e-5; %
    % [RegI]=CalcRegImg(bodycoil,nWrapI,wrap_m,Par);
    RegI = bodycoil*Par.Areg;
%     [RecMtx,g_map] = CalcReconMtrx(SSM,RegI,wrap_m,Ksi);
    [RecMtx,g_map] = CalcReconMtrx_sense_caipi(SSM,RegI,wrap_m,Ksi,Par,flag_CAIPI);
                
    if Par.nPE>nrPE
        [RecMtx] = Interp_RecMtx(RecMtx,Par.nPE,nFE,Par.nSL);
    end
end
%-----------------------------------------------------------------
%reconstruction
%-----------------------------------------------------------------
% 
[CompI]=ReconSENSE(RecMtx,nWrapI);


Inew = CompI + PriI;
%abs(CompI).*exp(i*Par.phase) + PriI;
