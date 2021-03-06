function [RegI]=CalcRegImg(PriI,WrapI,wrap_m,Par)
%This function is used to calculate the regularization image
%Please NOTICE RegI = 1./R
%INPUT
%   Pri_I: prior information of the reconstructed image [nPE nFE]
%   wrap_m: the matrix show overlaping positions [nPE/Rf nFE]
%   Par.Areg: weight for regularization
%   Par.alpha: regularization power of upper part
%   Par.beta: regularization power of lower part
%   Par.Kar: regularization threshold
%OUTPUT
%   RegI: Regularization image, [nPE nFE]
%Reference
%   RS:Reconstruction Page 232~234
%   Notation is same as the reference
%Note
%   Created by Feng on Aug 12th, 2008
PriI = abs(PriI); %only absolute value is used for regularization
[nPE nFE] = size(PriI);
%-------------------
%   Calculate Rm
%-------------------
indHigh = find(PriI > Par.Kar);
indLow = find(PriI <= Par.Kar);
Rm = zeros(nPE,nFE);
Rm(indHigh) = PriI(indHigh).^Par.alpha;
Rm(indLow) = PriI(indLow).^Par.beta.*(Par.Kar^(Par.alpha-Par.beta));
%-------------------
%   Calculate components of Rf
%-------------------
%Calculate reduection factor matrix
Sf = size(wrap_m,2)/size(wrap_m,1);
%Calculate SAm
SAm = WrapI.*conj(WrapI);
SAm = sum(SAm(:));
%Calculate SApc
SApc = PriI.*conj(PriI);
SApc = sum(SApc(:));
%-------------------
%   Calculate RegI
%-------------------

RegI = SApc*Sf/SAm./Rm;
RegI(isinf(RegI))=1e+5;
RegI = RegI*Par.Areg;
% RegI = Sf*Par.Areg./Rm;