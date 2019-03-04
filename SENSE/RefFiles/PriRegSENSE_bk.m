function [Inew]=PriRegSENSE(Kdata,PriI,Isens,Ksi,Par)
%This function is for the regularized mSENSE
%This version is only for 2D images
%INPUT
%   Kdata: partial k-space data [nPE nFE nCH]
%   PriI: Prior information of the image, a 3D matrix
%   Isens: sensitivity maps [nPE nFE nCH]
%   Ksi : noise correlation matrix
%   Par: parameters
%       .priLmd: weight Lamdba for the prior informaiton
%OUTPUT
%   Inew: reconstructed image

%----------------------------------------------------------
%       Preparation
%----------------------------------------------------------
[nPE,nFE,nCH]=size(Kdata);
[wrap_m] = sense_wrap_matrix(nPE,Par.Rf,1);
for ch=1:nCH
    PK = zeros(nPE,nFE);
    PK(Par.st_pt:Par.Rf:nPE,:) = Kdata(Par.st_pt:Par.Rf:nPE,:,ch);
    wrapped_image(:,:,ch) = K2Img(PK)*Par.Rf;
end
Lmd = Par.priLmd;
red_nrow = size(wrap_m,1);
invKsi = inv(Ksi);
%----------------------------------------------------------
%       Reconstruction
%----------------------------------------------------------
%******Set initial value*****
Inew=zeros(nPE,nFE);
% g_map = ones(nPE,nFE);
%******reconstruction part****
for row = 1:red_nrow
    ind_row = find(wrap_m(row,:));
    num_p = max(size(ind_row));
    for col =1:nFE
        AP = conj(squeeze(Isens(ind_row,col,:)));
        if num_p==1
            AP = AP.';
        end
        b = squeeze(wrapped_image(row,col,:));
        X0 = PriI(ind_row,col);  
        sumA = sum(abs(AP),2);
        indH = find(sumA>0);
        AP = AP(indH,:);
        X0 = X0(indH);
        X = zeros(num_p,1);
%         g_tmp = ones(num_p,1);
        %-------------------------------------
        %  Define matrices
        %-------------------------------------
        if ~isempty(AP)
            B = AP*invKsi;
            A_wave = B*AP';
            A_head = A_wave'*A_wave + Lmd^2*eye(max(size(indH)));%(num_p);%
            if rem(nPE,Par.Rf)>0
                inv_A_head = inv(A_head);
                ReconMtr = zeros(num_p,nCH);
                ReconMtr(indH,:) = inv_A_head*A_wave'*B;
                nb = squeeze(wrapped_image(ind_row,col,:));
                X = sum(ReconMtr.*nb,2);
                X(indH) = X(indH)+ inv_A_head*((Lmd^2)*X0);
            else
                b_wave = B*b;
                b_head = A_wave'*b_wave + Lmd^2*X0;
                X(indH) = inv(A_head)*b_head;
            end
        end
        %-------------------------------------
        %  Reconstruction
        %-------------------------------------       
        Inew(ind_row,col) = X;%inv(A_head)*b_head;%
%         g_tmp(indH)= sqrt(diag(inv(A_head'*A_head)).*diag(A_head'*A_head));
%         g_map(ind_row,col)=g_tmp;
    end %for col
end %for row


