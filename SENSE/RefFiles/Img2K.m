function [K]=Img2K(Img)
%This function is used to change image into k-space

% K = fftshift(fft2(Img)); %all philips data, GE head
% K = fft2(ifftshift(Img)); %%Shepp-Logan Phantom
% K = fftshift(fft2(ifftshift(Img,2))); 
% K = (fft2(ifftshift(Img,2)));
% K = fftshift(ifft2(fftshift(Img))); %cardiac, CAI, Shepp-Logan Phantom, Qu Peng
 K = fftshift(fft2(fftshift(Img))); %cardiac, CAI, Philips spine
% K = ifftshift(fft2(ifftshift(Img,1)));%each read out slice of 3D data