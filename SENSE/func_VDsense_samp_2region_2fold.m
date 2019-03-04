function [ParK_und_region1, ParK_und_region2, alias_m_region1, alias_m_region2] = func_VDsense_samp_2region_2fold(ParK_0,nPE_act,fac_region1,Rf_y,Rf_z,flag_opt,alias_m_region1,alias_m_region2)

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
Rf_y_region1 = Rf_y(1); %1.5;%1,1.5;%1.6,1.5;%2;2.4;3
Rf_z_region1 = Rf_z(1);%1;

if rem(Rf_y_region1,1)==0
    flag_opt = 1;
else
    flag_opt = 2;
end
% figure,imshow(abs(squeeze(ParK_0(:,12,:,1))),[]);imcontrast

switch flag_opt
    case 1
        ParK = ParK_0(1:Rf_y_region1:end,1:Rf_z_region1:end,:,:);
        [nrPE, nrSL,nFE, nCH] = size(ParK);
        ParK_v2_zpad = ParK;
        clear ParK;
    case 2
         %produce aliasing matrix
%         [alias_m]=sense_alias_matrix([nPE nSL],[nPE/Rf_y nSL/Rf_z]);
        if nargin<7
            [alias_m_region1]=sense_alias_matrix([nPE nSL],[nPE/Rf_y_region1 nSL/Rf_z_region1]);
        end
%         figure,imshow(alias_m_region1{1})
%         figure,imshow(alias_m_region1{2})
        %produce the full FOV wrapped image
%         [WrapI] = genWrapI(ParK_0,alias_m);            
        for ch=1:nCH
            for FE =1:nFE
                sWrapI(:,:,FE,ch) = alias_m_region1{1}*squeeze(tmp_img_full(:,:,FE,ch))*alias_m_region1{2};
            end
        end
%         as(sos(sWrapI,4))
        sWrapI_gen = permute(sWrapI,[1 3 2 4]);
        for ch=1:nCH
            ParK_v1(:,:,:,ch) = permute(Img2K_3d(permute(sWrapI(:,:,:,ch),[1 3 2 4])),[1 3 2 4]);
        end 
%         figure,imshow(abs(squeeze(ParK_v1(:,12,:,1))),[0 620]);%imcontrast
        ParK_v2 = crop(ParK_v1,[nPE_act/Rf_y_region1,nSL,nFE,nCH]);
%         figure,imshow(abs(squeeze(ParK_v2(:,12,:,1))),[]);imcontrast
        ParK_v2_zpad = zpad(ParK_v2,[nPE/Rf_y_region1,nSL,nFE,nCH]);
%         figure,imshow(abs(squeeze(ParK_v2_zpad(:,12,:,1))),[0 620]);%imcontrast
%         [nrPE, nrSL,nFE, nCH] = size(ParK_v2);
%         for ch=1:nCH
%             sWrapI_gen_2(:,:,:,ch) = K2Img_3d (permute(ParK_v2_zpad(:,:,:,ch),[1 3 2 4]));
%         end
%         as(sWrapI_gen_2)
end

ParK_und_region1 = zpad(crop(ParK_v2_zpad,[nPE_act/Rf_y_region1*fac_region1,nSL,nFE,nCH]),[nPE/Rf_y_region1,nSL,nFE,nCH]);
% figure,imshow(abs(squeeze(ParK_und_region1(:,12,:,1))),[0 620]);%imcontrast

%% Region-2
ParK_region2_fromRegion1 = ParK_v2_zpad - ParK_und_region1;
% figure,imshow(abs(squeeze(ParK_region2_fromRegion1(:,12,:,1))),[0 620]);%imcontrast

Rf_y_region2 = Rf_y(2)/Rf_y(1); 
Rf_z_region2 = Rf_z(2)/Rf_z(1);
ParK_und_region2 = ParK_region2_fromRegion1(1:Rf_y_region2:end,1:Rf_z_region2:end,:,:);
% figure,imshow(abs(squeeze(ParK_und_region2(:,12,:,1))),[0 620]);%imcontrast

end
