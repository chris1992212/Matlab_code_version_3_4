%% Generate mask from full reconstruction Img
% for i=1:11
%     folderpath = ['G:\WorkingFile\MiaoCodes\data\DATA', num2str(i),'\'];
%     load([folderpath, 'full_recon.mat'],'reconImg' );
%     mask=CreateMask3d(abs(squeeze(reconImg(:,:,:,1))));
%     save([folderpath, 'full_recon.mat'],'mask','-append' );
%     disp(folderpath);
% end
%% check mask
% for i=1:11
%     folderpath = ['G:\WorkingFile\MiaoCodes\data\DATA', num2str(i),'\'];
%     load([folderpath, 'full_recon.mat'],'mask' );
%     as(mask);
% end