%% Preprocess the raw data
% clear all
% Data_Preprocessing('E:\Qin_recon\data\', 4, 11 )
%% Reconstruction of the simulated PF+SENSE data
for i=2:3
    folderpath = ['E:\Qin_recon\data\DATA', num2str(i),'\']
    TestNum_s = 18;
    TestNum_e = 34;
    PFSENSE_main(folderpath, TestNum_s, TestNum_e);
end