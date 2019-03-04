function  Data_Preprocessing(filepath, start_num, end_num )
%DATA_PREPROCESSING is used to pre-process the original and save new data in the disk.
addpath('files\ImageDenoising\')
addpath('files\SENSE\')
addpath('files\SensitivityMap\')
addpath('files\parF\'); % new files created by Miao
addpath('files\utils\'); % new files created by Miao
addpath('files\ImageEvaluation\');% new files created by Miao
addpath('files\PFSENSE\');% new files created by Miao
if ~isdir(filepath)
    error('The first argument must be directory!');
    return
end

for i=start_num:end_num
    tempfile = strcat(filepath,'DATA',num2str(i), '\');
    load([tempfile ,'original.mat']);
    [raw3d_CoilImg,sense_map3d,body_coil3d, ~] = ...
    Multi_channel_address(imgDataCopy_FA1, imgDataCopy_FA2,PCICCMap_mpr, QBC_mpr);
    save(strcat(tempfile,'preprocessed.mat'),'raw3d_CoilImg','sense_map3d','body_coil3d','-v7.3');
    disp(['DATA ', num2str(i), ' has been saved!']);
    clear imgDataCopy_FA1 imgDataCopy_FA2 PCICCMap_mpr QBC_mpr

    disp('Start to recontruct the full FOV data...');
    acceleRatio.hf_fac_ky = 1;
    acceleRatio.hf_fac_kz = 1;
    acceleRatio.SENSEy_f = 1;
    acceleRatio.SENSEz_f = 1;
    acceleRatio.Regu = 1e-6;
    % load([folderpath,'preprocessed.mat'],'raw3d_CoilImg');
    PS_Img  = CreatePFSENSE(raw3d_CoilImg,acceleRatio); % return PFSENSE Img
    clear raw3d_CoilImg
    % load([folderpath,'preprocessed.mat'],'sense_map3d','body_coil3d');
    [FullReconImg, ~] = PFSENSERecon( PS_Img, sense_map3d, body_coil3d, acceleRatio);
    mask=CreateMask3d(abs(squeeze(FullReconImg(:,:,:,1))));
    save([tempfile,'full_recon.mat'],'FullReconImg','mask','-v7.3');
    savepath = [tempfile,'result\'];
    if ~exist(savepath,'dir')
        mkdir(savepath);
    end
end
