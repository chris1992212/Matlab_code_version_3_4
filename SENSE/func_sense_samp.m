function [ParK,alias_m,sWrapI_gen_2] = func_sense_samp(ParK_0,nPE_act,Rf_y,Rf_z,flag_opt,alias_m)

% %% Downsampling
% flag_opt = 2;
% Rf_y = 1.8;%1.6,1.5;%2;
% Rf_z = 1;
% nrPE = nPE/Rf_y;
% nrSL = nSL/Rf_z;

[nPE, nSL, nFE, nCH] = size(ParK_0);
switch flag_opt
    case 1
        ParK = ParK_0(1:Rf_y:end,1:Rf_z:end,:,:);
        % ParK = ParK_0(1:2:end,1:2:end,:,:);
        % ParK = ParK(1:4:end,1:2:end,:,:);
        [nrPE, nrSL,nFE, nCH] = size(ParK);
%         for ch=1:nCH
%             sWrapI_direct(:,:,:,ch) = K2Img_3d(permute(ParK(:,:,:,ch),[1 3 2 4]));
%         end
    case 2
        %produce aliasing matrix
        if nargin<6
            [alias_m]=sense_alias_matrix([nPE nSL],[nPE/Rf_y nSL/Rf_z]);
        end
%         figure,imshow(alias_m{1})
%         figure,imshow(alias_m{2})
        %produce the full FOV wrapped image
%         [WrapI] = genWrapI(ParK_0,alias_m);    
        for ch=1:nCH
            tmp_img_v2(:,:,:,ch) = K2Img_3d(permute(ParK_0(:,:,:,ch),[1 3 2 4]));
%             tmp_img_v2(:,:,:,ch) = K2Img_3d(permute(ParK_0_hf(:,:,:,ch),[1 3 2 4]));
        end 
        tmp_img_full = permute(tmp_img_v2,[1 3 2 4]);
        for ch=1:nCH
            for FE =1:nFE
                sWrapI(:,:,FE,ch) = alias_m{1}*squeeze(tmp_img_full(:,:,FE,ch))*alias_m{2};
            end
        end
        sWrapI_gen = permute(sWrapI,[1 3 2 4]);
        for ch=1:nCH
            ParK_v1(:,:,:,ch) = permute(Img2K_3d(permute(sWrapI(:,:,:,ch),[1 3 2 4])),[1 3 2 4]);
        end 
        ParK = zpad(crop(ParK_v1,[nPE_act/Rf_y,nSL,nFE,nCH]),[nPE/Rf_y,nSL,nFE,nCH]);  %% add on 2017/10/30
%         ParK = ParK_v1;
        [nrPE, nrSL,nFE, nCH] = size(ParK);
        for ch=1:nCH
            sWrapI_gen_2(:,:,:,ch) = K2Img_3d(permute(ParK(:,:,:,ch),[1 3 2 4]));
        end
%         as(sWrapI_gen_2)
end
% max(abs(sWrapI_direct(:))-abs(sWrapI_gen(:)))
% max(abs(ParK(:))-abs(ParK_v2(:)))
% figure,imshow(abs(abs(squeeze(ParK_v2(:,end/2,:,1)))-abs(squeeze(ParK(:,end/2,:,1)))),[]);colorbar;
% max(abs(sWrapI_direct(:))-abs(sWrapI_gen_2(:)))

% test_img = tmp_img_v2(:,:,24,1);
% figure,imshow(abs(test_img),[]);colorbar;
% test_img2K = Img2K_3d(test_img);
% test_img2K2img = K2Img_3d(test_img2K);
% figure,imshow(abs(test_img2K2img),[]);colorbar;
end

