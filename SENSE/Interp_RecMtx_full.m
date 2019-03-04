function [nRecMtx] = Interp_RecMtxhigh(RecMtx,nPE,nFE,nSL)
%This function is used to interpolation
if (nargin==3)
    %define target locations
    [XI,YI] = meshgrid(1:nFE,1:nPE);
    %define source grids
    [rnPE rnFE nCH] = size(RecMtx);
%     dPE = (nPE - 1)/rnPE;  
%     dFE = (nFE - 1)/rnFE;  
%     st_PE = dPE/2+1;
%     st_FE = dFE/2+1;
    
    dPE = (nPE + 6)/rnPE;
    dFE = (nFE + 6)/rnFE;
    st_PE = dPE/2-2;
    st_FE = dFE/2-2;
    
    [X,Y] = meshgrid(st_FE:dFE:(nFE + 6 - st_FE-dFE),st_PE:dPE: (nPE + 6 - st_PE - dPE));
    
    for ch=1:nCH
        [RealP] = interp2(X,Y,real(RecMtx(:,:,ch)),XI,YI);
        [ImgP] = interp2(X,Y,imag(RecMtx(:,:,ch)),XI,YI);
        nRecMtx(:,:,ch) = complex(RealP,ImgP);
    end
    nRecMtx(isnan(nRecMtx)) = 0;
    
elseif (nargin == 4)
    %define target locations
%     [XI,YI,ZI] = meshgrid(1:nSL,1:nFE,1:nPE);
 [XI,YI,ZI] = meshgrid(1:nSL,1:nPE,1:nFE); %%wyr change the order
    %define source locations
    [rnPE, rnSL, rnFE, nCH] = size(RecMtx);
    
    dPE = (nPE + 6)/rnPE;
    dFE = (nFE + 6)/rnFE;
    dSL = (nSL + 2)/rnSL;
    st_PE = dPE/2-2;
    st_FE = dFE/2-2;
    st_SL = dSL/2-1;
 
%     [X,Y,Z] = meshgrid(st_SL:dSL:nSL,st_FE:dFE:nFE,st_PE:dPE:nPE);
   [X,Y,Z] = meshgrid(st_SL:dSL: (nSL + 2 + st_SL - dSL),st_PE:dPE: (nPE + 6 + st_PE - dPE),st_FE:dFE:(nFE + 6 + st_FE-dFE)); %%wyr change the order
    
%     dPE = (nPE - 1)/rnPE;
%     dFE = (nFE - 1)/rnFE;
%     dSL = (nSL + 3)/rnSL;
%     st_PE = dPE/2+1;
%     st_FE = dFE/2+1;
%     st_SL = dSL/2-1;
%     [X,Y,Z] = meshgrid(st_SL:dSL: (nSL + 3 + st_SL - dSL),st_PE:dPE:nPE,st_FE:dFE:nFE); %%wyr change the order
   for ch=1:nCH
        [RealP] = interp3(X,Y,Z,real(RecMtx(:,:,:,ch)),XI,YI,ZI);
        [ImgP] = interp3(X,Y,Z,imag(RecMtx(:,:,:,ch)),XI,YI,ZI);
        nRecMtx(:,:,:,ch) = complex(RealP,ImgP);
    end
    nRecMtx(isnan(nRecMtx)) = 0;
end