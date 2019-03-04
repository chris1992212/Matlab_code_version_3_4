function [Img]= K2Img(K)

%  Img = ifft2(fftshift(K)); %all philips data, GE head
% Img = fftshift(ifft2(K)); %
% Img = fftshift(ifft2(fftshift(K)),2);
% Img = fftshift(ifft2((K)),2);
% Img = fftshift(fft2(fftshift(K))); %cAI, caridac, Shepp-Logan Phantom,Qu Peng
Img = ifftshift(ifft2(ifftshift(K))); %cAI, caridac,Philips spine
% Img = fftshift(ifft2(fftshift(K)),1);%3D data