function [senseMap,bodycoil,filled_mask]=SM_LresImg(Img, Par,bodycoil)
%this function is used to calculate sensitivity maps from low resolution
%images
%INPUTs
%   Img: low resolution images [nPE nFE nCH]
%   Par.flagSM: flag for smooth
%       0. without smoothing
%       1. use B-spline for smoothing
%OUTPUTs
%   senseMap:sensitivity maps in the same size as Img
%   bodycoi: the bodycoil image
%   filled_mask: mask of the image region
%----------------------------------------------------
%    Feng Huang, Aug 17, 2010
%----------------------------------------------------
ndim = size(Img);
if ndim ==3
    [nPE nFE nCH] = size(Img);
    % if bodycoi is not given, then SSoS is used
    if(nargin==2)
        %calculate bodycoil
        bodycoil = sqrt(sum(abs(Img).^2,3));
    end
    
    %find background
    [ind_big,Im]=sep_1d(abs(bodycoil(:)),1);
    mask = zeros(nPE,nFE);
    mask(ind_big{1})=1;
    se = strel('disk',2); %default is 2
    filled_mask = imfill(medfilt2(imdilate(double(mask),se)),'holes');
    
    %Calculate sensitivity maps
    indMask = find(filled_mask);
    senseMap = zeros(nPE,nFE,nCH);
    for ch =1:nCH
        tmpSM = zeros(nPE,nFE);
        tmpImg = Img(:,:,ch);
        tmpSM(indMask) = tmpImg(indMask)./bodycoil(indMask);
        senseMap(:,:,ch) = tmpSM;
    end
    
    %if required, then interpolate and extrapolate the sensitivity maps
    
    if Par.flagSM ==1
        [input_x_grid,input_y_grid] = meshgrid(floor(linspace(1,nFE,25)),...
            floor(linspace(1,nPE,25)));
        %only choose the reliable points
        ind_input = sub2ind([nPE nFE],input_y_grid(:),input_x_grid(:));
        ind_input = ind_input(mask(ind_input)>0);
        [input_x, input_y] = ind2sub([nPE nFE],ind_input);
        [x,y] = meshgrid(1:nFE,1:nPE);
        
        for ch=1:nCH
            raw_sensitivity = senseMap(:,:,ch);
            %smooth real part
            spline_function = tpaps([input_x(:),input_y(:)]',real(raw_sensitivity(ind_input))');
            val_spline = reshape(fnval(spline_function,[input_y_grid(:),input_x_grid(:)]'),size(input_x_grid));
            sensitivity_real = interp2(input_x_grid,input_y_grid,val_spline,x,y);
            %smooth imaginary part
            spline_function = tpaps([input_x(:),input_y(:)]',imag(raw_sensitivity(ind_input))');
            val_spline = reshape(fnval(spline_function,[input_y_grid(:),input_x_grid(:)]'),size(input_x_grid));
            sensitivity_imag = interp2(input_x_grid,input_y_grid,val_spline,x,y);
            %combine
            tmpSM = filled_mask.*complex(sensitivity_real,sensitivity_imag);
            senseMap(:,:,ch) = raw_sensitivity.*mask  + tmpSM.*(1-mask);
        end
    end
    
    %nomalize sensitivity maps
    SS = sqrt(sum(abs(senseMap).^2,3));
    
    for ch = 1:nCH
        tmpSM = zeros(nPE,nFE);
        tmpImg = senseMap(:,:,ch);
        tmpSM(indMask) = tmpImg(indMask)./SS(indMask);
        senseMap(:,:,ch) = tmpSM;
    end
else
    [nPE nFE nSL nCH] = size(Img);
    % if bodycoi is not given, then SSoS is used
    if(nargin==2)
        %calculate bodycoil
        bodycoil = sqrt(sum(abs(Img).^2,4));
    end
    
    %find background
    [ind_big,Im]=sep_1d(abs(bodycoil(:)),1);
    mask = zeros(nPE,nFE,nSL);
    mask(ind_big{1})=1;
    se = strel('disk',2); %default is 2
    filled_mask = imfill(imdilate(double(mask),se),'holes');
    
    %Calculate sensitivity maps
    indMask = find(filled_mask);
    senseMap = zeros(nPE,nFE,nSL,nCH);
    for ch =1:nCH
        tmpSM = zeros(nPE,nFE,nSL);
        tmpImg = Img(:,:,:,ch);
        tmpSM(indMask) = tmpImg(indMask)./bodycoil(indMask);
        senseMap(:,:,:,ch) = tmpSM;
    end
    
 
    if(nargin==2)
        %nomalize sensitivity maps
        SS = sqrt(sum(abs(senseMap).^2,4));
        
        for ch = 1:nCH
            tmpSM = zeros(nPE,nFE,nSL);
            tmpImg = senseMap(:,:,:,ch);
            tmpSM(indMask) = tmpImg(indMask)./SS(indMask);
            senseMap(:,:,:,ch) = tmpSM;
        end
    end
end



