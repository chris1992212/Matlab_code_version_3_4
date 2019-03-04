function [nRecMtx] = Interp_RecMtx(RecMtx,nPE,nFE,nSL)
%This function is used to interpolation
if (nargin==3)
    %define target locations
    [XI,YI] = meshgrid(1:nFE,1:nPE);
    %define source grids
    [rnPE rnFE nCH] = size(RecMtx);
    dPE = (nPE - 1)/rnPE;
    dFE = (nFE - 1)/rnFE;
    st_PE = dPE/2+1;
    st_FE = dFE/2+1;
    
    
    [X,Y] = meshgrid(st_FE:dFE:nFE,st_PE:dPE:nPE);
    
    for ch=1:nCH
        [RealP] = interp2(X,Y,real(RecMtx(:,:,ch)),XI,YI);
        [ImgP] = interp2(X,Y,imag(RecMtx(:,:,ch)),XI,YI);
        nRecMtx(:,:,ch) = complex(RealP,ImgP);
    end
    nRecMtx(isnan(nRecMtx)) = 0;
    
elseif (nargin == 4)
    %define target locations
%     [XI,YI,ZI] = meshgrid(1:nSL,1:nFE,1:nPE);
    [XI,YI,ZI] = meshgrid(1:nSL,1:nPE,1:nFE);
    %define source locations
    [rnPE, rnSL, rnFE, nCH] = size(RecMtx);
    
    dPE = (nPE - 1)/rnPE;
    dFE = (nFE - 1)/rnFE;
    dSL = (nSL - 1)/rnSL;
    st_PE = dPE/2+1;
    st_FE = dFE/2+1;
    st_SL = dSL/2+1;
    
    
%     [X,Y,Z] = meshgrid(st_SL:dSL:nSL,st_FE:dFE:nFE,st_PE:dPE:nPE);
    [X,Y,Z] = meshgrid(st_SL:dSL:nSL,st_PE:dPE:nPE,st_FE:dFE:nFE);
    
    for ch=1:nCH
        [RealP] = interp3(X,Y,Z,real(RecMtx(:,:,:,ch)),XI,YI,ZI);
        [ImgP] = interp3(X,Y,Z,imag(RecMtx(:,:,:,ch)),XI,YI,ZI);
        nRecMtx(:,:,:,ch) = complex(RealP,ImgP);
    end
    nRecMtx(isnan(nRecMtx)) = 0;
end