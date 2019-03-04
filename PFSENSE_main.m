 function PFSENSE_main(folderpath, TestNum_s, TestNum_e)
%This is the main script of the PF with SENSE algo for BrainQuant.
% It replaces the orginal SENSE_preprocessing function. 
% Created by Miao on 25 FEB 2019.
addpath('files\ImageDenoising\')
addpath('files\SENSE\')
addpath('files\SensitivityMap\')
addpath('files\parF\'); % new files created by Miao
addpath('files\utils\'); % new files created by Miao
addpath('files\ImageEvaluation\');% new files created by Miao
addpath('files\PFSENSE\');% new files created by Miao
% 
% folderpath = 'G:\WorkingFile\MiaoCodes\data\DATA1\';

%%%%%%% Generate self-defined Accelerating rate(Matrix) for PF and SENSE
acceleRatio = GenerateAccelerate();
%%%%%%% Load pre-process data
load([folderpath,'preprocessed.mat']);
%%%%% compared with the full reconstruction result
load([folderpath,'full_recon.mat']);
%%%%% Start Reconstruction mission
for i = TestNum_s:TestNum_e
    %%%%% Generate under-sampling PF+SENSE Img
    PS_Img  = CreatePFSENSE(raw3d_CoilImg,acceleRatio(i)); % return PFSENSE Img
    tic;
    [reconImg, G_map] = PFSENSERecon( PS_Img, sense_map3d, body_coil3d, acceleRatio(i));
    [ qc ] = QuantCriteria( reconImg, FullReconImg, logical(mask) );
    savepath = [folderpath,'result\'];
    save([savepath,'Test', num2str(acceleRatio(i).TestNum) , '.mat'],'reconImg','G_map','qc','-v7.3');
    toc;
end