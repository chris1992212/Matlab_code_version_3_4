function [Img]= K2Img_00(K)

% if ndims(K)==2
%     Img = ifftshift(ifft2(ifftshift(K))); %cAI, caridac,Philips spine
% %     Img = fftshift(ifft2(fftshift(K)),1);%3D data
% else
%     Img = ifftshift(ifftn(ifftshift(K))); %cAI, caridac,Philips spine
%     Img = fftshift(Img,2);
% end
 
S = size(K);
fctr = S(1)*S(2)*S(3);

K = reshape(K,S(1),S(2),S(3),prod(S(4:end)));

res = zeros(size(K));
for n=1:size(K,4)
%     res(:,:,:,n) = sqrt(fctr)*ifftn(ifftshift(x(:,:,:,n)));
%     res(:,:,:,n) = sqrt(fctr)*ifftshift(ifftn(fftshift(x(:,:,:,n))));
    res(:,:,:,n) = ifftshift(ifftn(fftshift(K(:,:,:,n))));
end

Img = reshape(res,S);
