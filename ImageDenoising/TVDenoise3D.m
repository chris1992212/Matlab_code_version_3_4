function [I,Out] = TVDenoise3D(varargin)
% This function  solves the TVL1-L2 model in 3D:
%
%   min aTV*TV(I) +  0.5|I - nI|_2^2
%INPUTS

% nI: the noisy image, a matrix
% aTV: weight for the TV term
% TVtype: 1 or 2, corresponding to anisotropic/isotropic TV in (2), recommend 2;
% opts      --- contains parameters for algorithm
%        * opts.maxItr: maximal total iteration number
%        * opts.gamma: 1.0 ~ 1.618
%        * opts.beta:  1 ~ 100
%        * opts.relchg_tol: stopping tolerance of norm(U-U_previous,'fro')/norm(U,'fro')
%        * opts.normalize: whether or not normalizes parameters and data
% 	     * opts.weight: the spatial weight for image denoising

% OUTPUTS: 
%          
%  I      -- denoised image
%  Out    -- iteration information, e.g., iter number, relative errors,
%            function values, etc.
% opts.weight = abs(g_map-1)*mean(abs(g_map(:)));
%---------------------------------------------------------------------
% Feng Huang, Aug 16, 2010
% This is a modification of RecRF by Yin Zhang and Wotao Yin from Rice Univ
% Invivo Corporation, Copyright (2010)
%---------------------------------------------------------------------


%%       Define Parameters
[nI,opts,URange] = ParseInputs(varargin{:});
aTV = opts.aTV;
TVtype = opts.TVtype;
maxItr = opts.maxItr;
gamma = opts.gamma;
beta = opts.beta;
relchg_tol = opts.relchg_tol;
Wts = opts.weight;

[nPE, nFE, nSL] = size(nI);
bComplex = true;    % use complex computation or not

%% normalize parameters and data
if opts.normalize
    if ~isreal(URange)||~isscalar(URange)||URange<eps; error('URange must be postive and real.'); end
    fctr = 1/URange;
    nI = fctr*nI;

    %aTV = sqrt(nPE*nFE)*aTV;
end

%% initialize constant parts of numinator and denominator (in order to save computation time)
I = nI;
B = fftn(nI);
Numer1 = B; %sqrt(nPE*nFE)*

Denom1 = ones(nPE,nFE,nSL); 
%change convolution in image space to be multiplication in k-space
prd = sqrt(aTV*beta);
Dz(1,1,1) = prd;
Dz(1,1,2) = -prd;
Denom2 = abs(psf2otf([prd,-prd],[nPE,nFE, nSL])).^2 + ...
    abs(psf2otf([prd;-prd],[nPE,nFE, nSL])).^2+...
    abs(psf2otf(Dz,[nPE,nFE, nSL])).^2;

Denom = Denom1 + Denom2;

%% initialize constants
[Ux,Uy,Uz]=Compute_Ux_Uy_Uz(nI);
% bx = complex(zeros(nPE,nFE));
% by = complex(zeros(nPE,nFE));
if isreal(nI)
    bx = zeros(nPE,nFE,nSL); 
    by = zeros(nPE,nFE,nSL);
    bz = zeros(nPE,nFE,nSL);
else
    bx = complex(zeros(nPE,nFE,nSL));
    by = complex(zeros(nPE,nFE,nSL));
    bz = complex(zeros(nPE,nFE,nSL));
end
%%  Main loop
for ii = 1:maxItr
        
    % ================================
    %  Begin Alternating Minimization
    % ----------------
    %   W-subprolem
    % ----------------
    switch TVtype
        case 1;   % anisotropic TV
            Ux = Ux + bx; 
            Uy = Uy + by;      % latest Ux and Uy are already calculated
            Uz = Uz + bz;
            Wx = sign(Ux).* max(abs(Ux)-Wts/beta,0);
            Wy = sign(Uy).* max(abs(Uy)-Wts/beta,0);
            Wz = sign(Uz).* max(abs(Uz)-Wts/beta,0);
        case 2;   % isotropic TV
            [Wx, Wy, Wz] = Compute_Wx_Wy_Wz_PtWs(Ux,Uy,Uz,bx,by,bz,Wts/beta);
        otherwise; 
            error('TVtype must be 1 or 2');
    end

    

    % ----------------
    %   U-subprolem
    % ----------------
    Iprev = I;
    
    rhs = Compute_rhs_DxtU_DytU_DztU(Wx,Wy,Wz,bx,by,bz,(aTV*beta)); 
    

    I = ifftn((Numer1 + fftn(rhs))./Denom); 
    
    if ~bComplex; I=real(I); end
    
    [Ux,Uy,Uz]=Compute_Ux_Uy_Uz(nI);
    
    %
    %  End Alternating Minimization
    % ================================

    % -------------------------------------------------
    % check stopping criterion
    %
    relchg = norm(I(:)-Iprev(:),'fro')/norm(I(:),'fro');

    if relchg < relchg_tol
        break;
    end
    
    % ------------------------------------------
    % Bregman update
    %
    bx = bx + gamma*(Ux - Wx);
    by = by + gamma*(Uy - Wy);
    bz = bz + gamma*(Uz - Wz);

end % outer

Out.iter = ii;
%% reverse normalization
if opts.normalize; I = I/fctr; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parse inputs
%%%

function [nI,opts,URange] = ParseInputs(varargin)

% error(nargchk(1,2,nargin,'struct'));
switch nargin
    case 1       % TVDenoise(nI)
        nI = varargin{1};
        opts.aTV = 0.01;
        opts.TVtype = 2;
        opts.maxItr = 100;       % max # of iterations
        opts.gamma = 1.0;        % noisy choice = 1.0
        opts.beta = 10;          % noisy choice = 10
        opts.relchg_tol = 5e-4;  % stopping tolerance based on relative change
        [nPE,nFE,nSL] = size(nI);
        opts.weight = ones(nPE,nFE,nSL);
        opts.normalize = 1;
    case 2       % TVDenoise(nI,opts)
        nI = varargin{1};
        opts = varargin{2};
        if ~isfield(opts, 'aTV')
            opts.aTV = 0.01; 
        end
        if ~isfield(opts, 'TVtype')
            opts.TVtype = 2; 
        end
        if ~isfield(opts, 'maxItr')
            opts.maxItr = 100; 
        end
        if ~isfield(opts, 'gamma')
            opts.gamma = 1; 
        end
        if ~isfield(opts, 'beta')
            opts.beta = 10; 
        end
        if ~isfield(opts, 'relchg_tol')
            opts.relchg_tol = 5e-4; 
        end
        if ~isfield(opts, 'weight')
            [nPE,nFE,nSL] = size(nI);
            opts.weight = ones(nPE,nFE,nSL);
        end
        if ~isfield(opts, 'normalize')
            opts.normalize = 1;
        end
end;

URange = max(abs(nI(:)));

