function [U, Out] = LSBB(V, B, A, opts)
% min aP/2*||U-V||^2 + 1/2*||A*U-B||^2
%  U
% Linearize the second term using BB step
% Xiaojing Ye, Spring 2010.

%% Parameters
aP = opts.aP;
innitr = opts.innitr;
eTol = opts.eTol;

%% Initialization
iter = 0; U = V; OBJ = realmax; RelErr = 1; CV = aP*V; bb = 1; N=numel(V);

%% Inner loop for LS
while RelErr > eTol && iter < innitr

    Err = A*U - B; BB = norm(Err(:))^2/N;
    Grad = A'*Err;

    if iter > 0
        DiffU = U - Uprev;
        DiffGrad = Grad - Gradprev;
        Numer = AtimesB(DiffGrad,DiffU);
        bb = Numer/norm(DiffU,'fro')^2;
%         fprintf('BB Step Size: %.3f\n',bb)
        if (bb < 1e-10) || (bb > 1e+10)
            fprintf('Oops, BB stepsize so small/big: %e\n',bb)
            bb = 2;    
        end
    end;

    OBJprev = OBJ; 
%     OBJ1 = norm(U-V,'fro')^2; OBJ2=BB;
    OBJ = .5*(aP*norm(U-V,'fro')^2 + BB);
    RelErr = abs(OBJ-OBJprev)/OBJ;
    %fprintf('InnItr=%d, OBJ1=%.3f, OBJ2=%.3f, OBJ=%.3f\n',iter, OBJ1, OBJ2, OBJ);
    Uprev = U;

    % Update U
    U = (bb*U + CV - Grad)/(bb+aP);
%     fprintf('Unorm=%.3f, Vnorm=%.3f, Grad=%.3f\n',max(abs(U(:))),max(abs(V(:))),max(abs(Grad(:))));

    iter = iter + 1;
    
    Gradprev = Grad;   
end
Out.BB = BB;
return;