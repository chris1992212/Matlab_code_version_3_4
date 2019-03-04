function [I,Out] = TVDenoise(varargin)
% This function  solves the TVL1-L2 model:
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
% Feng Huang, Aug 11, 2010
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

[m, n] = size(nI);
bComplex = true;    % use complex computation or not

%% normalize parameters and data
if opts.normalize
    if ~isreal(URange)||~isscalar(URange)||URange<eps; error('URange must be postive and real.'); end
    fctr = 1/URange;
    nI = fctr*nI;

    %aTV = sqrt(m*n)*aTV;
end

%% initialize constant parts of numinator and denominator (in order to save computation time)
I = nI;
B = fft2(nI);
Numer1 = B; %sqrt(m*n)*
Denom1 = ones(m,n); 
%change convolution in image space to be multiplication in k-space
prd = sqrt(aTV*beta);
Denom2 = abs(psf2otf([prd,-prd],[m,n])).^2 + abs(psf2otf([prd;-prd],[m,n])).^2;

Denom = Denom1 + Denom2;

%% initialize constants
[Ux,Uy]=Compute_Ux_Uy(nI);
% bx = complex(zeros(m,n));
% by = complex(zeros(m,n));
if isreal(nI)
    bx = zeros(m,n); by = zeros(m,n);
else
    bx = complex(zeros(m,n));
    by = complex(zeros(m,n));
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
            Ux = Ux + bx; Uy = Uy + by;      % latest Ux and Uy are already calculated
            Wx = sign(Ux).* max(abs(Ux)-Wts/beta,0);
            Wy = sign(Uy).* max(abs(Uy)-Wts/beta,0);
        case 2;   % isotropic TV
            [Wx, Wy] = Compute_Wx_Wy_PtWs(Ux,Uy,bx,by,Wts/beta);
            %[Wx, Wy] = Compute_Wx_Wy(Ux,Uy,bx,by,1/beta);
        otherwise; 
            error('TVtype must be 1 or 2');
    end

    

    % ----------------
    %   U-subprolem
    % ----------------
    Iprev = I;
    
    rhs = Compute_rhs_DxtU_DytU(Wx,Wy,bx,by,(aTV*beta)); 
    

    I = ifft2((Numer1 + fft2(rhs))./Denom); 
    
    if ~bComplex; I=real(I); end
    
    [Ux,Uy]=Compute_Ux_Uy(I);
    
    %
    %  End Alternating Minimization
    % ================================

    % -------------------------------------------------
    % check stopping criterion
    %
    relchg = norm(I-Iprev,'fro')/norm(I,'fro');

    if relchg < relchg_tol
        break;
    end
    
    % ------------------------------------------
    % Bregman update
    %
    bx = bx + gamma*(Ux - Wx);
    by = by + gamma*(Uy - Wy);

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
        [m,n] = size(nI);
        opts.weight = ones(m,n);
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
            [m,n] = size(nI);
            opts.weight = ones(m,n);
        end
        if ~isfield(opts, 'normalize')
            opts.normalize = 1;
        end
end;

URange = max(abs(nI(:)));

