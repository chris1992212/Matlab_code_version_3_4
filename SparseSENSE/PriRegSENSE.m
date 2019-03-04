function [Inew]=PriRegSENSE(Kdata,PriI,Isens,Ksi,Par)
%This function is for the prior information regularized SENSE
%This version is for both 2D and 3D images
%INPUT
%   Kdata: partial k-space data [rnPE nFE nCH]
%   PriI: Prior information of the image, a 2D matrix, [nPE nFE]
%   Isens: sensitivity maps [nPE nFE nCH]
%   Ksi : noise correlation matrix
%   Par: parameters
%       .priLmd: weight Lamdba for the prior informaiton
%OUTPUT
%   Inew: reconstructed image
%-----------------------------------------------------------------
%               Feng Huang, Aug 13, 2010
%      Based on Eqs. [3] and [6] in Fa-Hsuan Lin's paper in 2004
%   1. One midification is let image be 0 at locations all sensitivity maps
%   are 0 to reduce g-factor
%   2. To reduce reconstruction time, g-map is not calculated
%Aug 24, 2010: The 3D version is added
%-----------------------------------------------------------------

if ndims(Kdata)==3 %2D SENSE
    
    %----------------------------------------------------------
    %       Preparation
    %----------------------------------------------------------
    
    % Prepare some parameters
    
    [red_nrow,nFE,nCH]=size(Kdata);
    
    nPE = Par.nPE;
    %produce aliasing matrix
    [wrap_m]=sense_alias_matrix(Par.nPE,red_nrow);
    %produce the reduced FOV wrapped image
    wrapped_image = zeros(red_nrow,nFE,nCH);
    for ch =1:nCH
        wrapped_image(:,:,ch) = K2Img(Kdata(:,:,ch));
    end
    
    Lmd = Par.priLmd;
    % prepare for whitening, this is for general Ksi. If Ksi is diagonal, this
    % step can simplified
    [Vn,D] = eig(Ksi);
    invKsi = diag(1./(sqrt(diag(D))));
    %----------------------------------------------------------
    %       Reconstruction
    %----------------------------------------------------------
    %******Set initial value*****
    Inew=zeros(nPE,nFE);
    %******reconstruction part****
    for row = 1:red_nrow
        ind_row = find(wrap_m(row,:));
        num_p = max(size(ind_row));
        for col =1:nFE
            AP = conj(squeeze(Isens(ind_row,col,:)));
            A = AP';
            
            A = invKsi*Vn'*A;
            
            if num_p==1
                A = A.';
            end
            b = squeeze(wrapped_image(row,col,:));
            b = invKsi*Vn'*b;
            X0 = PriI(ind_row,col);
            sumA = sum(abs(A),1);
            indH = find(sumA>0);
            A = A(:,indH);
            X0 = X0(indH);
            X = zeros(num_p,1);
            %-------------------------------------
            %  Define matrices
            %-------------------------------------
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
            %-------------------------------------
            %  Reconstruction
            %-------------------------------------
            Inew(ind_row,col) = X;%inv(A_head)*b_head;%
        end %for col
    end %for row
elseif ndims(Kdata)==4 %3D SENSE
     %----------------------------------------------------------
    %       Preparation
    %----------------------------------------------------------
    
    % Prepare some parameters
    
    [rnPE,rnSL,nFE,nCH]=size(Kdata);
    
    nPE = Par.nPE;
    nSL = Par.nSL;
    %produce aliasing matrix
    [wrap_m]=sense_alias_matrix([nPE nSL],[rnPE rnSL]);
    %produce the reduced FOV wrapped image
    wrapped_image = zeros(rnPE,rnSL,nFE,nCH);
    for ch =1:nCH
        wrapped_image(:,:,:,ch) = K2Img(Kdata(:,:,:,ch));
    end
    
    Lmd = Par.priLmd;
    % prepare for whitening, this is for general Ksi. If Ksi is diagonal, this
    % step can simplified
    [Vn,D] = eig(Ksi);
    invKsi = diag(1./(sqrt(diag(D))));
    %----------------------------------------------------------
    %       Reconstruction
    %----------------------------------------------------------
    %******Set initial value*****
    Inew=zeros(nPE,nSL,nFE);
    %******reconstruction part****
    for SL = 1:rnSL
        ind_sl = find(wrap_m{2}(:,SL));
        num_sl = max(size(ind_sl));
        for PE = 1:rnPE
            ind_row = find(wrap_m{1}(PE,:));
            num_p = max(size(ind_row));
            Rf = num_p*num_sl;
            for FE =1:nFE
                AP = conj(squeeze(Isens(ind_row,ind_sl,FE,:)));
                AP = reshape(AP,Rf,nCH);
                A = AP';
                
                A = invKsi*Vn'*A;
                
%                 if Rf==1
%                     A = A.';
%                 end
                b = squeeze(wrapped_image(PE,SL,FE,:));
                b = invKsi*Vn'*b;
                X0 = PriI(ind_row,ind_sl,FE);
                X0 = reshape(X0,Rf,1);
                sumA = sum(abs(A),1);
                indH = find(sumA>0);
                A = A(:,indH);
                X0 = X0(indH);
                X = zeros(Rf,1);
                %-------------------------------------
                %  Define matrices
                %-------------------------------------
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
                %-------------------------------------
                %  Reconstruction
                %-------------------------------------
                Inew(ind_row,ind_sl,FE) = reshape(X,num_p,num_sl);%inv(A_head)*b_head;%
            end %for FE
        end %for PE
    end %for SL
end

