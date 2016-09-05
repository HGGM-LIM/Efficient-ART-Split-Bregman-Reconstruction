function u   = ARTReconstruction(A,d,relaxParam,numIterART,u)
% u   = ARTReconstruction(A,d,relaxParam,numIterART,u)
%
% Algebraic reconstruction technique, also known as Kaczmarz method
%
% Inputs:
% 
% A             = System matrix or Jacobian matrix, nr x nc
%               = d, data, nr x 1
% relaxParam    = relaxation parameter (between 0 and 1). Choose 0.1 to
% remove noise and 0.9 to get a closer fit of the data
% numIterART    = number of iterations
% u             = initial guess of the solution image, nc x 1. Initialize
% as zeros(nc,1) 
%
% Outputs:
%
% u             = solution image, nc x 1
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

% Norm of rows
Anorm               = sum(A.*A,2);

% Loop across iterations
for it = 1:numIterART
    % Loop across rows
    for ir = 1:size(A,1)
        num         = d(ir)-A(ir,:)*u;
        u           = u + (relaxParam*num/Anorm(ir))*A(ir,:)';
    end % ir
end % it