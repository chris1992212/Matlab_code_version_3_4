function [cv, C] = ncov(I,ind)
%this function is to calculate noise correlation matrix
%INPUTS
%   I: a 3D matrix [nPE nFE nCH]
%   ind: the index of location of noise 
%OUTPUTS
%   cv: noise covariance matrix
%   C:  noise correlation matrix
nmtx = [];
for ch=1:size(I,3)
    tmpI = I(:,:,ch);
    nmtx = [nmtx tmpI(ind)];
end

cv=cov(double(nmtx));
if (nargout > 1)
 C=corrcoef(double(nmtx));
end
