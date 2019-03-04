% this script is used to merge the preprocessing of training data into one function
% input: 
% Matname: Struct of volunteer data
% y_factor: SENSE factor in phase direction
% Areg:regulazation 
function [PS_Img, Par] = SENSE_Preprocessing(raw3d_CoilImg,acceleRatio)

PS_Img  = CreatePFSENSE(raw3d_CoilImg,acceleRatio); % return PFSENSE Img
end

