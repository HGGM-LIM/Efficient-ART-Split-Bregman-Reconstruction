function u = TV_SB_denoising_2D(f,mu,lambda,nin,nout)
% u = TV_SB_denoising_2D(f,mu,lambda,alpha,nin,nout)
%
% Solves the 2D-TV denoising problem using the SplitBregman formulation
%
% Inputs:
%
% f         = noisy image, 2D, nx x ny
% temporal dimensions)
% mu        = 1, weight of data constraint. Increase for faster
% convergence, decrease for noisy data (1 usually works fine)
% lambda    = 1, weight of TV constraint
% alpha     = 1, sparsity parameter weighting TV functional
% nInner    = 5 to 10. Number of inner iterations that impose TV
% constraints. The higher the number of iterations the more is imposed the
% TV constraint
% nOuter    = 2. Number of outer iterations that impose data fidelity term.
% Choosing this larger than 1 can correct for a decrease in contrast that
% TV may provide
%
% Split Bregman parameters (mu, lambda, alpha) allow to tune the algorithm
% if needed but it is generally quite robust and selecting all parameters
% to 1 (or as suggested in the paper) works fine. Then, the only parameters
% that need to be chosen are the number of iterations.
%
% Outputs:
%
% u         = denoised image
%
% The method is computationally efficient for large scale problems, as it
% requires only solving a linear system, which is solved in the Fourier
% domain without increasing the size of the operators (which remain of size
% nx xn ny), and solving thresholding operations. This makes it one of
% the fastest and computationally efficient state-of-the art methods.
%
%
% This code is a modified version of spatial TV reconstruction Goldstein'n
% code mrics.m downloaded from
% (http://www.ece.rice.edu/~tag7/Tom_Goldstein/Split_Bregman.html), see Tom
% Goldstein and Stanley Osher. The Split Bregman Method for L1-Regularized
% Problems. SIAM J. Imaging Sci., 2(2), 323–343.
%
%
% Code downloaded from the repository
% https://github.com/HGGM-LIM/Efficient-ART-Split-Bregman-Reconstruction
%
% If you use this code, please cite Chamorro-Servent et al. Use of Split
% Bregman denoising for iterative reconstruction in fluorescence diffuse
% optical tomography. J Biomed Opt, 18(7):076016, 2013.
% http://dx.doi.org/10.1117/1.JBO.18.7.076016
%
% Juan FPJ Abascal, Judit Chamorro-Servent, Juan Aguirre
% Departamento de Bioingenieria e Ingenieria Aeroespacial
% Universidad Carlos III de Madrid, Madrid, Spain
% juanabascal78@gmail.com, juchamser@gmail.com, desco@hggm.es

[rows,cols] = size(f);

% Image normalization to be independent of the selection of mu. It is
% multiplied by 10 so there is output of the shrinkage on the first
% iterations (faster convergence)
normFactor  = 10/max(abs(f(:)));
f           = normFactor*f;

% Reserve memory for the auxillary variables
u           = f;
f0          = f;
x           = zeros(rows,cols);
y           = zeros(rows,cols);
bx          = zeros(rows,cols);
by          = zeros(rows,cols);

% Hessian operator in the Fourier domain
uker        = zeros(rows,cols);
uker(1,1) = 4; uker(1,2) = -1; uker(2,1) = -1;
uker(rows,1) = -1; uker(1,cols) = -1;
uker        = fft2(uker);
uker        = mu + lambda*uker;

% Outer loop
for outer = 1:nout
    for inner=1:nin
        % Inner loop
        rhs     = mu*f+lambda*Dxt(x-bx)+lambda*Dyt(y-by);
        
        % Solve the linear system in the Fourier domain
        u       = ifft2(fft2(rhs)./uker);
        
        % Compute gradients of the image
        dx      = Dx(u);
        dy      = Dy(u);
        
        x       = shrink1(dx+bx, 1/lambda);
        y       = shrink1(dy+by, 1/lambda);
        %[x,y] = shrink2( dx+bx, dy+by,1/lambda);
        
        % Bregman parameters for TV contraint
        bx      = bx+dx-x;
        by      = by+dy-y;
    end
    % Bregman parameter for data fidelity term
    f = f + f0-u;
    
end

u   = u/normFactor;

return


function d = Dx(u)
[rows,cols] = size(u);
d = zeros(rows,cols);
d(:,2:cols) = u(:,2:cols)-u(:,1:cols-1);
d(:,1) = u(:,1)-u(:,cols);
return

function d = Dxt(u)
[rows,cols] = size(u);
d = zeros(rows,cols);
d(:,1:cols-1) = u(:,1:cols-1)-u(:,2:cols);
d(:,cols) = u(:,cols)-u(:,1);
return

function d = Dy(u)
[rows,cols] = size(u);
d = zeros(rows,cols);
d(2:rows,:) = u(2:rows,:)-u(1:rows-1,:);
d(1,:) = u(1,:)-u(rows,:);
return

function d = Dyt(u)
[rows,cols] = size(u);
d = zeros(rows,cols);
d(1:rows-1,:) = u(1:rows-1,:)-u(2:rows,:);
d(rows,:) = u(rows,:)-u(1,:);
return


function [xs,ys] = shrink2(x,y,lambda)

s = sqrt(x.*conj(x)+y.*conj(y));
ss = s-lambda;
ss = ss.*(ss>0);

s = s+(s<lambda);
ss = ss./s;

xs = ss.*x;
ys = ss.*y;

return;


function xs = shrink1(x,lambda)

s = abs(x);
xs = sign(x).*max(s-lambda,0);

return;
