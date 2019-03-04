function [Par] = default_par(Par)
%This function is used to define default parameters

%---------------------------------------------
% parameters for acquired data
%---------------------------------------------
if ~isfield(Par,'Rf') %reduction factor
    disp('It is necessary to know the reduction factor. Please Input!')
    exit; %should wait for the input and continue
end
if ~isfield(Par,'st_pt') 
    disp('It is assumed that the starting point of PE is 1')
    Par.st_pt = 1; %starting PE line
end

%---------------------------------------------
% parameters for CS
%---------------------------------------------
if ~isfield(Par,'aTV') 
    Par.aTV = 1e-2; %parameters for TV term
end
if ~isfield(Par,'aL1') 
    Par.aL1 = 0;%parameter for wavelet term
end
if ~isfield(Par,'Lg') 
    Par.Lg = 1.2; %no change at locations where g-factor is smaller than this value
end
%---------------------------------------------
% parameters for prior informatoin regularized SENSE
%---------------------------------------------
if ~isfield(Par,'priLmd')
    Par.priLmd = 0.2; %weight for regularization
end
if ~isfield(Par,'innitr') 
    Par.innitr = 15; %number of iternal iterations
end
if ~isfield(Par,'eTol')
    Par.eTol = 1e-4; %determine if it is converged
end