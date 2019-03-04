function [Inew,RecMtx]=PriRegSENSEtmp(Kdata,PriI,Isens,SSM,Ksi,Par)
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
%-----------------------------------------------------------------

%-----------------------------------------------------------------
%Produce wrapped image with reduced image support reduction
%-----------------------------------------------------------------
[nacqPE,nacqSL,nFE,nCH]=size(Kdata);

nPE = Par.nPE;
nSL = Par.nSL;

%------ Produce the full FOV wrapped image------ 

%produce aliasing matrix
[alias_m]=sense_alias_matrix([Par.nPE Par.nSL],[nacqPE nacqSL]);
%produce the full FOV wrapped image
[WrapI] = genWrapI(Kdata,alias_m);


ImgC = zeros(size(Isens));
for ch=1:nCH
    ImgC(:,:,:,ch) = PriI.*Isens(:,:,:,ch);
end
[nWrapI] = genWrapI(ImgC,alias_m);

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
%     [dx,dy] = gradient(abs(PriI));
    bodycoil = ones(nPE,nSL,nFE);%./(1e-2+Par.weight); %sqrt(dx.^2+dy.^2)
%     bodycoil = Par.weight;%./sqrt(abs(dx).^2+abs(dy).^2);
%     bodycoil =1./sqrt(abs(dx).^2+abs(dy).^2);
%     bodycoil(isnan(bodycoil))=1e+5;

    %acquired number of PE for lower resolution
    [nrPE,nrSL,nrFE,nCH]=size(SSM);
    
     %acquired number of PE for lower resolution
    ParNrPE = round(nacqPE/Par.nPE*nrPE);
    ParNrSL = round(nacqSL/Par.nSL*nrSL);
    [wrap_m]=sense_alias_matrix([nrPE nrSL],[ParNrPE ParNrSL]);
    
    [bodycoil] = abs(Interp_RecMtx(bodycoil,nrPE,nrFE,nrSL));
    
    Par.Areg = 200;% 1e-3;%5e-5; %
    % [RegI]=CalcRegImg(bodycoil,nWrapI,wrap_m,Par);
    RegI = bodycoil*Par.Areg;
    [RecMtx,g_map] = CalcReconMtrx(SSM,RegI,wrap_m,Ksi);
    
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
