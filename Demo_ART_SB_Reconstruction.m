% Demo_ART_SB_Reconstruction.m
% 
% Demo for a novel iterative reconstruction method that alternates a
% computationally efficient linear solver (ART) with a fast denoising step
% based on the Split Bregman formulation and that is used to reconstruct
% fluorescence diffuse optical tomography data, as in the paper J
% Chamorro-Servent, J F P J Abascal, J Aguirre, S Arridge, T Correia, J
% Ripoll, M Desco, J J Vaquero. Use of Split Bregman denoising for
% iterative reconstruction in fluorescence diffuse optical tomography. J
% Biomed Opt, 18(7):076016, 2013. http://dx.doi.org/10.1117/1.JBO.18.7.076016       
% 
% Code downloaded from the repository 
% https://github.com/HGGM-LIM/Efficient-ART-Split-Bregman-Reconstruction
%
% This work proposes to combine computationally efficient solver and
% denoising steps with the potential to handle large scale problems in a
% fast an efficient manner. Here we combine Algebraic Reconstruction
% Technique (ART) and Split Bregman denoising but these methods can be
% substituted by the method of choice. 
%
% Our particular choices are explained as follows: 
%    ART has been widely used in tomography to handle large data sets, as
% it can work with one data point (projection) at a time. 
% 
%    The Split Bregman formulation allows to solve the total variation
% denoising problem in a fast and computationally efficient way. At each
% iteration the solution is given by analytical formulas. A linear system  
% can be solved in the Fourier domain, keeping the size of the problem (the
% Hessian) equal to the number of voxels in the image, or by using
% Gauss-Seidel method, which exploits the block diagonal structure of the
% Hessian. 
%
% In this demo we solve the image reconstruction problem of fuorescence
% diffuse optical tomography but it can be applied to other imaging
% modalities. 
%
% 
% If you use this code, please cite Chamorro-Servent et al. Use of Split
% Bregman denoising for iterative reconstruction in fluorescence diffuse
% optical tomography. J Biomed Opt, 18(7):076016, 2013.
% http://dx.doi.org/10.1117/1.JBO.18.7.076016        
%
% Judit Chamorro-Servent, Juan FPJ Abascal, Juan Aguirre
% Departamento de Bioingenieria e Ingenieria Aeroespacial
% Universidad Carlos III de Madrid, Madrid, Spain
% juchamser@gmail.com, juanabascal78@gmail.com, desco@hggm.es

 
% READ SIMULATED DATA 
% Load data, Jacobian matrix and target image
load('DataRed','data','JacMatrix','uTarget');
 
[nr nc]     = size(JacMatrix);

% Display target image in the discretized domain
N           = [20 20 10];   
x           = linspace(-10,10,N(1));
y           = linspace(-10,10,N(2));
z           = linspace(0,10,N(3));
[X,Y,Z]     = meshgrid(x,y,z);
% Uncomment below to get a plot of all slices
% h           = Plot2DMapsGridSolution(reshape(uTarget,N),X,Y,Z,3); colorbar;

rand('seed',0);

% Add Gaussian noise
data       = data + 0.01*max(data)*randn(length(data),1);

% Randomized index of Jacobian rows for ART reconstruction
% indRand     = randperm(length(data));
% indRand     = 1:nr; % Uncomment to remove randomized index

% -------------------------------------------------------------------------
% ART
%
% It requires selecting few parameters: 
%    numIterART = number of iterations. The more iterations the better fit
% of the data (in this case it converges in 30 iterations) 
%    relaxParam = relaxation parameter (between 0 and 1). Choose 0.1 to
% remove noise and 0.9 to get a closer fit of the data
numIterART  = 30;   % 
relaxParam  = 0.1;
u0          = zeros(nc,1);

fprintf('Solving ART reconstruction ... (it takes around 1 s)\n');
tic; 
uART = ARTReconstruction(JacMatrix,data,relaxParam,numIterART,u0); 
toc
uART     = reshape(uART,N);
%  h = Plot2DMapsGridSolution(uART,X,Y,Z,3); colorbar;
       
% -------------------------------------------------------------------------
% ART-SB  
%
% It requres selecting the following paramters for ART and Split Bregman 
% denoising:   
%    ART (see previous section). Here we provide high relaxation for a
% better fit of the data, as noise will be removed in the denoising step) 
%
%    SPLIT BREGMAN DENOISING:
% Split Bregman parameters (mu, lambda, alpha) allow to tune the algorithm
% if needed but it is generally quite robust and selecting all parameters
% to 1 (or as suggested in the paper) works fine. Then, the only parameters
% that need to be chosen are the number of iterations 
%
%    mu     = weight of the fidelity term. (Generally 0.1 to 1 works fine;
%    see the paper for more details). The larger the more weight to the
%    noisy image
%    lamdba = weight of the functionals that impose TV constraints. The
% larger the higher the regularization. (Usually chosen larger than mu and
% any value around 1 works fine) 
%    alpha  = the weight of the functional that imposes TV sparsity
% (L1-norm). No need to tune this parameter, value 1 should work
%    nInner = Number of inner iterations that impose TV constraints. The
% higher the number of iterations the more is imposed the TV constraint
%    nOuter = Number of outer iterations that impose data fidelity
% constraint. Choose this larger than 1 (2 works fine), as TV may provide
% an output image with lower constrast and the Bregman iteration can
% correct for that. 
numIter     = 30;
numIterART  = 10;      
relaxParam  = 0.9;
uARTSB      = zeros(N);
mu          = 0.3;
lambda      = 2*mu;
alpha       = 1;
nInner      = 5;
nOuter      = 2; 
% indRand = randperm(length(data));
fprintf('Solving ART-SB reconstruction ... (it takes around 12 s)\n');
h = waitbar(0,'Solving ART-SB reconstruction') ;
tic
for it = 1:numIter
    % ART reconstruction step: Iterative linear solver
    sol = ARTReconstruction(JacMatrix,data,relaxParam,numIterART,uARTSB(:)); 
    solGrid     = reshape(sol,N);   

    % SB denoising step
    uARTSB      = TV_SB_denoising_3D(solGrid,mu,lambda,alpha,nInner,nOuter);    
    % Uncomment below to do 2D slice-by-slice smoothing instead. It takes
    % similar time but it could be parallelized, which can be faster in some
    % applications for large scale problems (now SB denoising takes less
    % than a second), by changing the for to a parfor loop  
%     for iz = 1:N(3)
%         uARTSB(:,:,iz) = TV_SB_denoising_2D(solGrid(:,:,iz),mu,lambda,nInner,nOuter);
%     end

    % Compute solution error norm
    err(it) = norm(uARTSB(:)-uTarget(:))/norm(uTarget(:));
    
    waitbar(it/numIter);
end % it
toc
close(h);

% Display results
figure; plot(err); ylabel('Solution error'); xlabel('Number of iterations');
title('Convergence of ART-SB');

% Target image
Plot2DMapsGridSolution(reshape(uTarget,N),X,Y,Z,3); 
set(gcf,'name','TARGET','numbertitle','off'); 
colormap gray;

% Reconstructed images0
Plot2DMapsGridSolution(uART,X,Y,Z,3); 
set(gcf,'name','ART reconstruction','numbertitle','off') 
colormap gray;

Plot2DMapsGridSolution(uARTSB,X,Y,Z,3); 
set(gcf,'name','ART-SB reconstruction','numbertitle','off') 
colormap gray;  
% -------------------------------------------------------------------------    
%
