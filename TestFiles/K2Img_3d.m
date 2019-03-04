function [Img]= K2Img_3d(K)

if ndims(K)==2
    Img = ifftshift(ifft2(ifftshift(K))); %cAI, caridac,Philips spine
%     Img = fftshift(ifft2(fftshift(K)),1);%3D data
else
    Img = ifftshift(ifftn(ifftshift(K))); %cAI, caridac,Philips spine
%     Img = fftshift(Img,2);
end