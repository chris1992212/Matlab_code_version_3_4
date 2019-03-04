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
    PK = Kdata(1:Par.Rf:nPE,:,ch);
    wrapped_image(:,:,ch) = K2Img(PK);
end
Lmd = Par.priLmd;
red_nrow = size(wrapped_image,1);
%----------------------------------------------------------
%       Reconstruction
%----------------------------------------------------------
%******Set initial value*****
Inew=zeros(nPE,nFE);
%******reconstruction part****
for row = 1:red_nrow
    ind_row = find(sign(wrap_m(row,:))==1);
    num_p = max(size(ind_row));
    for col =1:nFE
        A = zeros(nCH,num_p);
        b = zeros(nCH,1);
        rowA = 0;
        for ch =1:nCH
            SM = Isens(ind_row,col,ch);
            rowA = rowA + 1;
            A(rowA,:) = SM.';
            b(rowA) = wrapped_image(row,col,ch);             
        end %for ch
        X0 = PriI(ind_row,col);  
        %-------------------------------------
        %  deal with zero sensitivity
        %-------------------------------------
        sumA = sum(abs(A),1);
        indH = find(sumA>0);
        A = A(:,indH);
        X0 = X0(indH);
        
        X = zeros(num_p,1);
        
        if ~isempty(A)
            [U,S,V] = svd(A);   
            f = diag(S).^2./(diag(S).^2+ Lmd^2);
            sz = numel(indH);
            X_c = zeros(sz,1);  
            
            for j=1:sz
                X_c = X_c + (f(j)*U(:,j)'*b/S(j,j)+(1-f(j))*V(:,j)'*X0)*V(:,j);
            end %for j
            X(indH) = X_c;
        end
        Inew(ind_row,col) = X;
    end %for col
end %for row


