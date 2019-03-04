function  [MSE,PSNR] = MSE(image,denoiseimg, mask)
if nargin ==2
    mask = ones(size(image));
end
maximum = double(max(image(:)));
image = double(image) / maximum * 255;
denoiseimg = double(denoiseimg) / maximum * 255;
if(size(image) ~= size(denoiseimg))
    error('the size of these two inputs must be the same');
else

    %%%%%%%%%%
    MSE = sum(sum((mask.*(image-denoiseimg)).^2))/ sum(mask(:));
    PSNR = 10*log10(255^2 / MSE);
    %%%%%%%%%%
   
%   disp(['MSE  = ',num2str(MSE)]);
%   disp(['PSNR = ',num2str(PSNR)]);

end
end