function [K]= Img2K(Img)

if ndims(Img)==2
    K = fftshift(fft2(fftshift(Img))); %cAI, caridac,Philips spine
%     Img = fftshift(ifft2(fftshift(K)),1);%3D data
else
    K = fftshift(fftn(fftshift(Img))); %cAI, caridac,Philips spine
%     Img = fftshift(Img,2);
end