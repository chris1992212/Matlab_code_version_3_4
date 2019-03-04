function [WrapI] = genWrapI(ParK,wrap_m)
%This function is used to produce wrapped image from partially acquired
%data. It works for artitrary reduction factor
if ~iscell(wrap_m) %1D SENSE
    [nrPE nFE nCH]= size(ParK);
    sWrapI = zeros(nrPE,nFE,nCH);
    nPE = size(wrap_m,2);
    %check if the input is k-space or image space
    if nrPE == size(wrap_m,1) %k-space
        for ch=1:nCH
            %     sWrapI(:,:,ch) = fftshift(K2Img(ParK(:,:,ch)),1);
            sWrapI(:,:,ch) = K2Img(ParK(:,:,ch));
        end   
    else % image space
        if nrPE == nPE
            for ch=1:nCH
                sWrapI(:,:,ch) = wrap_m*ParK(:,:,ch);
            end
        end
    end
    
    WrapI = zeros(nPE,nFE,nCH);
    for PE = 1:nrPE
        wPE = find(wrap_m(PE,:));
        Rf = max(size(wPE));
        WrapI(wPE,:,:) = repmat(sWrapI(PE,:,:),[Rf 1 1]);
    end
else %2D SENSE
    [nrPE nrSL nFE nCH]= size(ParK);
    nPE = size(wrap_m{1},2);
    nSL = size(wrap_m{2},1);
    sWrapI = zeros(size(wrap_m{1},1),size(wrap_m{2},2),nFE,nCH);
    
    %check if the input is k-space or image space
    if (nrPE == size(wrap_m{1},1))&& (nrSL == size(wrap_m{2},2))%k-space
        for ch=1:nCH
            sWrapI(:,:,:,ch) = K2Img(ParK(:,:,:,ch));
        end   
    else % image space
        if (nrPE == nPE)&&(nrSL == nSL)
            for ch=1:nCH
                for FE =1:nFE
                    sWrapI(:,:,FE,ch) = wrap_m{1}*squeeze(ParK(:,:,FE,ch))*wrap_m{2};
                end
            end
        end
    end
  
    
    
    WrapI = zeros(nPE,nSL,nFE,nCH);
    for SL = 1:size(wrap_m{2},2)
        wSL = find(wrap_m{2}(:,SL));
        RfSL = max(size(wSL));
        for PE = 1:size(wrap_m{1},1)
            wPE = find(wrap_m{1}(PE,:));
            RfPE = max(size(wPE));
            WrapI(wPE,wSL,:,:) = repmat(sWrapI(PE,SL,:,:),[RfPE RfSL 1 1]);
        end
    end
    
end