function [ParK_und_region1, ParK_und_region2, alias_m_region1, alias_m_region2] = func_VDsense_samp_2region_overlap(ParK_0,nPE_act,fac_region1,Rf_y,Rf_z,flag_opt,overlap_lines,filter_opt,alias_m_region1,alias_m_region2)

% %% Downsampling
% flag_opt = 2;
% Rf_y = 1.8;%1.6,1.5;%2;
% Rf_z = 1;
% nrPE = nPE/Rf_y;
% nrSL = nSL/Rf_z;

[nPE, nSL, nFE, nCH] = size(ParK_0);

for ch=1:nCH
    tmp_img_v2(:,:,:,ch) = K2Img_3d (permute(ParK_0(:,:,:,ch),[1 3 2 4]));
end
tmp_img_full = permute(tmp_img_v2,[1 3 2 4]);

%% Region-1
region1_lines_ky = nPE_act*fac_region1;
region1_area = nPE/2-region1_lines_ky/2+1:nPE/2+region1_lines_ky/2;
% ParK_region1_full = ParK_0(region1_area,:,:,:);
% ParK_region1_full_zpad = zpad(ParK_region1_full,[nPE,nSL,nFE,nCH]);
% chk_img = ifft3c(permute(ParK_region1_full_zpad,[1 3 2 4]));
% figure,imshow(rot90(sos(squeeze(chk_img(:,:,end/2,:)))),[])

Rf_y_region1 = Rf_y(1); %1.5;%1,1.5;%1.6,1.5;%2;2.4;3
Rf_z_region1 = Rf_z(1);%1;

switch flag_opt
    case 1
%         figure,imshow(abs(squeeze(ParK_region1_full_zpad(:,end/2,:,1))),[]);imcontrast
        ParK_region1_full = ParK_0(region1_area,:,:,:);
        ParK_region1_full_zpad = zpad(ParK_region1_full,[nPE,nSL,nFE,nCH]);
        if filter_opt==1
            filter = Tukey_mod1114( region1_lines_ky,overlap_lines+1,0,0);
%             figure,plot(filter)
            filter_fin = repmat(zpad(filter,[nPE,1]),[1,nSL,nFE,nCH]);
            ParK_region1_full_zpad_filter = ParK_region1_full_zpad.*filter_fin;
%             figure,imshow(abs(squeeze(filter_fin(:,end/2,:,1))),[]);imcontrast
%             figure,imshow(abs(squeeze(ParK_region1_full_zpad_filter(:,end/2,:,1))),[]);imcontrast
            ParK = ParK_region1_full_zpad_filter(1:Rf_y_region1:end,1:Rf_z_region1:end,:,:);
         else
            ParK = ParK_region1_full_zpad(1:Rf_y_region1:end,1:Rf_z_region1:end,:,:);
        end
        [nrPE, nrSL,nFE, nCH] = size(ParK);
        ParK_v2_zpad = ParK;
        clear ParK;
%         for ch=1:nCH
%             sWrapI_direct(:,:,:,ch) = K2Img_3d (permute(ParK(:,:,:,ch),[1 3 2 4]));
%         end
    case 2
        %produce aliasing matrix
%         [alias_m]=sense_alias_matrix([nPE nSL],[nPE/Rf_y nSL/Rf_z]);
        if nargin<7
            [alias_m_region1]=sense_alias_matrix([nPE nSL],[nPE/Rf_y_region1 nSL/Rf_z_region1]);
        end
%         figure,imshow(alias_m{1})
%         figure,imshow(alias_m{2})
        %produce the full FOV wrapped image
%         [WrapI] = genWrapI(ParK_0,alias_m);            
        for ch=1:nCH
            for FE =1:nFE
                sWrapI(:,:,FE,ch) = alias_m_region1{1}*squeeze(tmp_img_full(:,:,FE,ch))*alias_m_region1{2};
            end
        end
%         as(sWrapI)
        sWrapI_gen = permute(sWrapI,[1 3 2 4]);
        for ch=1:nCH
            ParK_v1(:,:,:,ch) = permute(Img2K_3d(permute(sWrapI(:,:,:,ch),[1 3 2 4])),[1 3 2 4]);
        end 
        ParK_v2 = crop(ParK_v1,[region1_lines_ky/Rf_y_region1,nSL,nFE,nCH]);
        ParK_v2_zpad = zpad(ParK_v2,[nPE/Rf_y_region1,nSL,nFE,nCH]);
%         [nrPE, nrSL,nFE, nCH] = size(ParK_v2);
%         for ch=1:nCH
%             sWrapI_gen_2(:,:,:,ch) = K2Img_3d (permute(ParK_v2_zpad(:,:,:,ch),[1 3 2 4]));
%         end
%         as(sWrapI_gen_2)
end
ParK_und_region1 = ParK_v2_zpad;
% figure,imshow(abs(squeeze(ParK_und_region1(:,end/2,:,1))),[]);imcontrast

%% Region-2
% ParK_region2_full = ParK_0 - ParK_region1_full_zpad;
% ParK_0 - zpad(crop(ParK_0(region1_area,:,:,:),[region1_lines_ky,nSL,nFE,nCH]),[nPE,nSL,nFE,nCH]);
% max(abs(ParK_region2_full(:)-ParK_region2_full_0(:)))
% figure,imshow(abs(ParK_region2_full(:,:,end/2,1)),[]);imcontrast
% chk_img2 = ifft3c(permute(ParK_region2_full,[1 3 2 4]));
% figure,imshow(rot90(sos(squeeze(chk_img2(:,:,24,:)))),[])

% flag_opt = 2;
Rf_y_region2 = Rf_y(2);%1.5;%1,1.5;%1.6,1.5;%2;2.4;3
Rf_z_region2 = Rf_z(2);

clear ParK alias_m sWrapI sWrapI_gen ParK_v1 ParK_v2 ParK_v2_zpad sWrapI_gen_2;
switch flag_opt
    case 1
        if filter_opt==1
            ParK_region2_full = ParK_0 - ParK_region1_full_zpad_filter;
        else
            region1_area_nonoverlap = nPE/2-(region1_lines_ky/2-overlap_lines)+1:nPE/2+(region1_lines_ky/2-overlap_lines);
            ParK_region1_full_zpad_nonoverlap = ParK_region1_full_zpad(region1_area_nonoverlap,:,:,:);
            ParK_region1_full_zpad_nonoverlap = zpad(ParK_region1_full_zpad_nonoverlap,[nPE,nSL,nFE,nCH]);
%             figure,imshow(abs(squeeze(ParK_region1_full_zpad_nonoverlap(:,end/2,:,1))),[]);imcontrast
%             ParK_region2_full = ParK_0 - ParK_region1_full_zpad;
            ParK_region2_full = ParK_0 - ParK_region1_full_zpad_nonoverlap;
        end
%         figure,imshow(abs(squeeze(ParK_region2_full(:,end/2,:,1))),[]);imcontrast

        ParK = ParK_region2_full(1:Rf_y_region2:end,1:Rf_z_region2:end,:,:);
        [nrPE, nrSL,nFE, nCH] = size(ParK);
        ParK_v2_zpad = ParK;
        clear ParK;
%         for ch=1:nCH
%             sWrapI_direct(:,:,:,ch) = K2Img_3d (permute(ParK(:,:,:,ch),[1 3 2 4]));
%         end
    case 2
        %produce aliasing matrix
%         [alias_m]=sense_alias_matrix([nPE nSL],[nPE/Rf_y nSL/Rf_z]);
        if nargin<8
            [alias_m_region2]=sense_alias_matrix([nPE nSL],[nPE/Rf_y_region2 nSL/Rf_z_region2]);
        end
%         figure,imshow(alias_m{1})
%         figure,imshow(alias_m{2})
%         %produce the full FOV wrapped image
%         for ch=1:nCH
%             tmp_img_v2(:,:,:,ch) = K2Img_3d (permute(ParK_0(:,:,:,ch),[1 3 2 4]));
%         end 
%         tmp_img_full = permute(tmp_img_v2,[1 3 2 4]);
        for ch=1:nCH
            for FE =1:nFE
                sWrapI(:,:,FE,ch) = alias_m_region2{1}*squeeze(tmp_img_full(:,:,FE,ch))*alias_m_region2{2};
            end
        end
%         as(sWrapI)
        sWrapI_gen = permute(sWrapI,[1 3 2 4]);
        for ch=1:nCH
            ParK_v1(:,:,:,ch) = permute(Img2K_3d(permute(sWrapI(:,:,:,ch),[1 3 2 4])),[1 3 2 4]);
        end 
%         figure,imshow(abs(ParK_v1(:,:,end/2,1)),[]);imcontrast
%         ParK_v2 = ParK_v1 - zpad(crop(ParK_v1,[region1_lines_ky/Rf_y_region2,nSL,nFE,nCH]),[nPE/Rf_y_region2,nSL,nFE,nCH]);
        ParK_v2 = zpad(crop(ParK_v1,[nPE_act/Rf_y_region2,nSL,nFE,nCH]),[nPE/Rf_y_region2,nSL,nFE,nCH]) ...
            - zpad(crop(ParK_v1,[region1_lines_ky/Rf_y_region2,nSL,nFE,nCH]),[nPE/Rf_y_region2,nSL,nFE,nCH]); % modified on 10/31
        ParK_v2_zpad = ParK_v2;%zpad(ParK_v2,[nPE/Rf_y_region1,nSL,nFE,nCH]);
%         figure,imshow(abs(ParK_v2_zpad(:,:,end/2,1)),[]);imcontrast
%         [nrPE, nrSL,nFE, nCH] = size(ParK_v2);
%         for ch=1:nCH
%             sWrapI_gen_2(:,:,:,ch) = K2Img_3d (permute(ParK_v2_zpad(:,:,:,ch),[1 3 2 4]));
%         end
%         as(sWrapI_gen_2)
end
ParK_und_region2 = ParK_v2_zpad;
% figure,imshow(abs(squeeze(ParK_und_region2(:,end/2,:,1))),[]);imcontrast



