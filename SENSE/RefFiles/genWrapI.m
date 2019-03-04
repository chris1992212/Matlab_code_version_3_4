function [WrapI] = genWrapI(ParK,wrap_m)
%This function is used to produce wrapped image from partially acquired
%data. It works for artitrary reduction factor
[nrPE nFE nCH]= size(ParK);
nPE = size(wrap_m,2);
sWrapI = zeros(nrPE,nFE,nCH);
for ch=1:nCH
    sWrapI(:,:,ch) = fftshift(K2Img(ParK(:,:,ch)),1);
%     sWrapI(:,:,ch) = K2Img(ParK(:,:,ch));
end
WrapI = zeros(nPE,nFE,nCH);
for PE = 1:nrPE
    wPE = find(wrap_m(PE,:));
    Rf = max(size(wPE));
    WrapI(wPE,:,:) = repmat(sWrapI(PE,:,:),[Rf 1 1]);
end