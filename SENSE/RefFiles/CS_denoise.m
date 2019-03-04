function [nI]=CS_denoise(CompI,g_map,Par)
%This function is used to reduce noise level in CompI
%INPUT
%   CompI: the image reconstructed by SENSE
%   g_map: the g_factor map 
%OUTPUT
%   nI: cleared new image
%fftshift k-space for the fast algorithm
ImgK = fftshift(Img2K(CompI)); 
U0 = ifft2(ImgK);
gK = fftshift(Img2K(g_map));
g_factor = ifft2(gK);
g_factor = abs(g_factor)./min(abs(g_factor(:)));

[m, n] = size(ImgK);

B = ImgK(:);
W = [];
WT = [];
%------------------------------------------------
% Define parameters
%------------------------------------------------
aTV = Par.aTV;   % 1e-4 is good for noisy (sigma = .01) under-sampled phantom

%%---------------------
% added for wavelets
aL1 = Par.aL1;%1e-2;
if aL1>0
    addpath('rwt');
    wav = daubcqf(2);
    W = @(x) midwt(x,wav);
    WT = @(x) mdwt(x,wav);
end
%%---------------------

opts = [];
opts.idisp = 1;
g_factor(find(abs(g_factor)<1))=1;
opts.Ig = abs(g_factor);
mask = zeros(m,n);
mask(1:Par.Rf:m,:)=1;
mask = fftshift(mask);
opts.acq = find(mask==1);
clear g_factor gK

pick=true(m,n);

%------------------------------------------------
% Call reconstruction function
%------------------------------------------------
[nI,Out_RecPF] = RecPF(m,n,aTV,aL1,pick,B,2,opts,WT,W,range(abs(U0(:))),U0);
nI = fftshift(nI);