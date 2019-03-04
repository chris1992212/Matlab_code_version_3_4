function [Recon]=ReconSENSE(RecMtx,WrapI)
%This function is used reconstrut the last image
% WrapI = repmat(WrapI,[Par.Rf 1 1]); 
if ndims(WrapI)==3
    Recon = sum(RecMtx.*WrapI,3);
else
    Recon = sum(RecMtx.*WrapI,4);
end