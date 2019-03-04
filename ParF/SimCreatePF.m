function [ ParK ] = SimCreatePF( ParK,hf_fac_ky,hf_fac_kz )
%The function is called to generate the partial fourier K-data for further simulation.
% The default acceleration ratio is 0.8 in ky and 0.9 in kz.
% Originally created by Aiqi Sun. Modified by Zhenzhuang Miao on 11 Feb 2019.

nPE_act = 244;%192;%144; % The phase encoding of sub-coils is 244 instead of 324 in the latest system. Zero is filled in both sides.
[nPE, nSL, ~, ~, ~] = size(ParK); 
nSL_act = nSL;
acq_lines_ky = ceil(nPE_act*hf_fac_ky); % round(nPE_act*hf_fac); 
acq_lines_kz = ceil(nSL_act*hf_fac_kz); % round(nPE_act*hf_fac);
offset_lines_ky = acq_lines_ky - nPE_act/2;
offset_lines_kz = acq_lines_kz - nSL_act/2;
ParK(nPE/2+offset_lines_ky+1:(nPE/2+nPE_act/2),:,:,:,:) = 0;
ParK(:,nSL/2+offset_lines_kz+1:(nSL/2+nSL_act/2),:,:,:) = 0;
end

