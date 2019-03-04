%=========================================================================
function sensitivity = CalcCoilSensitivities(data,bodycoil,mask);
% ------------------------------------------------------------------------
% Calculate sensitivity maps
% ------------------------------------------------------------------------

% ------------------------------------------------------------------------
% Get matrix sizes
x_res       = size(data,1);
y_res       = size(data,2);
nr_coils    = size(data,3);


% se = strel('disk',2); %default is 2
% interp_mask = imfill(medfilt2(imdilate(double(mask),se)),'holes');
% interp_mask = imfill(double(mask),'holes');
% ------------------------------------------------------------------------
% Smooth and interpolate sensitivity maps using a smoothing thin-plate
% spline. The real and imaginary parts of the sensitivity are smoothed
% separately as the smoothing function can only take real inputs.
for i = 1:nr_coils

    sensitivity(:,:,i) = mask.* data(:,:,i)./bodycoil;
    sensitivity(isnan(sensitivity)) = 0;
    
    
end
%only use the smoothed version to fill holes 

%normalization
SS = sqrt(sum(abs(sensitivity).^2,3));
for ch = 1:nr_coils
    sensitivity(:,:,ch) = sensitivity(:,:,ch)./SS;
end
 sensitivity(isnan(sensitivity)) = 0;