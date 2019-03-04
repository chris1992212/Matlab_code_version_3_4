function [K]=Img2K_00(Img)
%This function is used to change image into k-space

% K = fftshift(fft2(Img)); %all philips data, GE head
% K = fft2(ifftshift(Img)); %%Shepp-Logan Phantom
% K = fftshift(fft2(ifftshift(Img,2))); 
% K = (fft2(ifftshift(Img,2)));
% K = fftshift(ifft2(fftshift(Img))); %cardiac, CAI, Shepp-Logan Phantom, Qu Peng
%  K = fftshift(fft2(fftshift(Img))); %cardiac, CAI, Philips spine
% K = ifftshift(fft2(ifftshift(Img,1)));%each read out slice of 3D data

S = size(Img);
fctr = S(1)*S(2)*S(3);

Img = reshape(Img,S(1),S(2),S(3),prod(S(4:end)));

res = zeros(size(Img));
for n=1:size(Img,4)
% 	res(:,:,:,n) = 1/sqrt(fctr)*fftshift(fftn(x(:,:,:,n)));
%     res(:,:,:,n) = 1/sqrt(fctr)*fftshift(fftn(ifftshift(x(:,:,:,n))));
    res(:,:,:,n) = fftshift(fftn(ifftshift(Img(:,:,:,n))));
end

K = reshape(res,S);
