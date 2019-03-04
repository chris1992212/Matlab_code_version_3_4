function [Par]=define_Par()
% Initialize Parameters
%---------------------------------------------
% parameters for SENSE
%---------------------------------------------
%This function is used to define default values of parameters
%   Par.Areg: weight for regularization
%   Par.alpha: regularization power of upper part
%   Par.beta: regularization power of lower part
%   Par.Kar: regularization threshold
%   Par.ACS: definition of ACS lines
%           2d case: [ACS_st ACS_end]
%           3d case: [PE_st PE_end
%                     sl_st sl_end];  
%   Par.ind_filter: index of filter
%       1. hanning filter
%       2. filter for high pass 1d
%       3. filter for high pass 2d
%       4. half hanning filter
%       5. 2d hanning filter
%   Par.cw: definition of filter for high pass [c w]
%   Par.indSM: index for sensitivity maps
%       1. do nothing, use image as sensitivity maps
%       2. use sum of squares without smoothing
%       3. use sum of squares as denominator
%       4. use pesudo body coil as denominator
%   Par.flagSM: flag for smooth
%       0. without smoothing
%       1. use B-spline for smoothing

%fixed parameters
Par.Areg = 1e-6;%0.0001;
Par.alpha = 2;
Par.beta = 1;
Par.Kar = 0;

Par.Rf = 4; %acceleration factor
Par.st_pt = 1; %starting PE line


%---------------------------------------------
% parameters for CS
%---------------------------------------------
Par.aTV = 1e-2; %parameters for TV term
Par.aL1 = 0;%parameter for wavelet term
Par.Lg = 1.2; %no change at locations where g-factor is smaller than this value
%---------------------------------------------
% parameters for prior informatoin regularized SENSE
%---------------------------------------------

Par.priLmd = 0.4; %weight for regularization
Par.innitr = 5; %number of iternal iterations
Par.eTol = 1e-4; %determine if it is converged