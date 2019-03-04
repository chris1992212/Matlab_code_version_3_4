function [ acceleRatio ] = GenerateAccelerate()
%This function is used for generating a set of accelerating ratio for PF and SENSE in two
%sides.
[num,txt,raw] =xlsread('TestList.xls');
acceleRatio = cell2struct(raw(2:end,:),txt,2);

% acceleRatio.hf_fac_ky = 1;%0.8
% acceleRatio.hf_fac_kz = 1;%0.9
% acceleRatio.SENSEy_f = 2;% 2 
% acceleRatio.SENSEz_f = 1.5;% 1

end

