function ParK_Img  = CreatePFSENSE( ParK_Img,acceleR)
%The function is used for creating a set of simulated undersampling k-data from full
%measurement. It includes partial fourier data and SENSE accelerated in PE and SL.
% return argument is Img.
disp('The current accelerating rate is:');
disp(['SENSE y = ' num2str(acceleR.SENSEy_f) ...
    '; SENSE z = ' num2str(acceleR.SENSEz_f)]);
disp(['ParF y = ' num2str(acceleR.hf_fac_ky) ...
    '; ParF z = ' num2str(acceleR.hf_fac_kz)]);

[nPE, nFE, nSL, nCH, nEcho] = size(ParK_Img);
ParK_Img = fft3c(crop(ParK_Img,[324 nFE nSL nCH nEcho])); % Img->KData
ParK_Img = circshift(ParK_Img,[-1 0 0 0 0]);
ParK_Img = permute(ParK_Img,[1 3 2 4 5]); % [324 384 48 16 6]
%%%%%%%%%%%%%%%%% Downsampling for PF simulation
ParK_Img = SimCreatePF( ParK_Img,acceleR.hf_fac_ky,acceleR.hf_fac_kz );
%%$%%%%%%% Downsampling for SENSE in two directions   
if acceleR.SENSEy_f >= 1 || acceleR.SENSEz_f >= 1
     ParK_Img  = SimCreateSENSE( ParK_Img, acceleR.SENSEy_f, acceleR.SENSEz_f); %return Img
else
    error('The acclerating rate of the SENSE must be greater than 1!');
end
disp('Simulated PF+SENSE image has been generated!');
end

