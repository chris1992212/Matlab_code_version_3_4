function [nRecMtx] = Interp_RecMtx(RecMtx,nPE,nFE)
%This function is used to interpolation
%define target locations
[XI,YI] = meshgrid(1:nFE,1:nPE);
%define source locations
[rnPE rnFE nCH] = size(RecMtx);
if nPE >= rnPE
%     dPE = round((nPE - 1)/rnPE);
%     dFE = round((nFE - 1)/rnFE);
%     
%     st_PE = rem(nPE-1,dPE)/2+1;
%     st_FE = rem(nFE-1,dFE)/2+1;
    
    dPE = (nPE - 1)/rnPE;
    dFE = (nFE - 1)/rnFE;
    st_PE = dPE/2+1;
    st_FE = dFE/2+1;
   
else
    dPE = (nPE - 1)/rnPE;
    dFE = (nFE - 1)/rnFE;
    st_PE = dPE/2+1;
    st_FE = dFE/2+1;
%     st_PE = round((rnPE - nPE)/2);
%     ed_PE = st_PE + nPE -1;
%     st_FE = round((rnFE - nFE)/2);
%     ed_FE = st_FE + nFE -1;
%     for ch =1:nCH
%         Ktmp = fftshift(fft2(RecMtx(:,:,ch)));
%         nRecMtx(:,:,ch) = ifft2(fftshift(Ktmp(st_PE:ed_PE,st_FE:ed_FE)));
%     end
end

[X,Y] = meshgrid(st_FE:dFE:nFE,st_PE:dPE:nPE);

for ch=1:nCH
    [RealP] = interp2(X,Y,real(RecMtx(:,:,ch)),XI,YI);
    [ImgP] = interp2(X,Y,imag(RecMtx(:,:,ch)),XI,YI);
    nRecMtx(:,:,ch) = complex(RealP,ImgP);
end
nRecMtx(isnan(nRecMtx)) = 0;