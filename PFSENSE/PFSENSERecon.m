function [ reconI_final, G_map ] = PFSENSERecon( PS_Img,sense_map3d,body_coil3d, acceleR)
%-------------------------------------------------------------------
% Define parameters
%-------------------------------------------------------------------
% Rf = 4;
% DPA_recon = 0;%
% Rf_y = SENSEy_f;%324/130; %FA1£∫324/130; FA2:3£ª (~244/98) %324/116;(~244/87,limit:R=2.78/0.8/0.9=3.8)  % 2(144/72), 144/60, 144/54£¨ 144/48=3
Par.necho = 6;
Par.nFE = 384;
Par.nPE = 324;
Par.nSL = 48;
%parameters for Philips SENSE
% Par.Rf = Rf;
Par.Areg = acceleR.Regu;%1e-6;
Par.alpha = 2;
Par.beta = 1;
Par.Kar = 0;
%parameters for sensitivity maps
Par.flagSM = 0;
%parameters for image denoising
Par.aTV = 0.01;%0.02;
Par.TVtype = 1;%2;
Par.maxItr = 100;       
Par.gamma = 1.0;        
Par.beta = 10;          
Par.relchg_tol = 5e-4;  
Par.normalize = 1;
%parameters for prior information regularized SENSE
Par.priLmd = 0.1; %weight for regularization
Par.Lg = 1.0; %no change at locations where g-factor is smaller than this value
%%%%%%%%%%%%%%%%%%SENSE reconstruction%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------
% Initial Reconstruction
%-------------------------------------------------------------------
%     close all;
%     Par.nPE = nPE;%162;%192;
%     Par.nSL = nSL;%36;
%     Par.Areg = areg;%10;%0.1;%10; 1e-1;%1e-3;%10;%1e-6;%120;%10,50,100,120,150-200(Ω”Ω¸∆Ê“Ï),1e-6;%0.0001;
if ~exist('Phi')
    Phi = eye(16);
end
disp(['Start to reconstruct SENSE for regu=',num2str(Par.Areg)]);

reconI_final = zeros(Par.nPE, Par.nSL, Par.nFE, Par.necho,'single');
% RecMtx = zeros(Par.nPE, Par.nSL, Par.nFE, Par.necho,'single');
G_map = zeros(Par.nPE, Par.nSL, Par.nFE, Par.necho,'single');
for echo = 1:Par.necho
    [CompI,recmtx,g_map]=RegPreSENSE(PS_Img(:,:,:,:,echo),Phi,sense_map3d,body_coil3d,Par);
    reconI_final(:,:,:,echo)=single(CompI);%nFE,nPE,nSL
%     RecMtx(:,:,:,echo)=single(recmtx);%nFE,nPE,nSL
    G_map(:,:,:,echo)=single(g_map);%nFE,nPE,nSL
    % g_map_new=permute(g_map,[3 1 2]);
    % recon_3dsense_Rx_reg10 = CompI; g_map_Rx=g_map;
end
reconI_final = permute(reconI_final,[3,1,2,4]);
% RecMtx = permute(RecMtx,[3,1,2,4]);
G_map = permute(G_map,[3,1,2,4]);
end

