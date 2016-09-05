# Efficient-ART-Split-Bregman-Reconstruction
MATLAB demo to use a novel and efficient iterative reconstruction method that alternates the Algebraic Reconstruction Technique (ART) with a Total Variation (TV) denoising step based on the Split Bregman formulation.
 
The method is **ideal for large scale problems** as it proposes to combine an efficient linear solver with and an efficient denoising method. Here we use ART and the Split Bregman formulation, but these methods can be substituted by the method of choice. 

Our particular choices are explained as follows:

- ART has been widely used in tomography to handle large data sets, as it can work with one data point (projection) at a time. 

- The Split Bregman formulation allows to solve the total variation denoising problem in a fast and computationally efficient way. At each iteration the solution is given by analytical formulas. A linear system  can be solved in the Fourier domain, keeping the size of the problem (the Hessian) equal to the number of voxels in the image, or by using Gauss-Seidel method, which exploits the block diagonal structure of the Hessian. 

This demo reconstructs fluorescence diffuse optical tomography data as presented in the paper **Use of Split Bregman denoising for iterative reconstruction in fluorescence diffuse optical tomography. J Chamorro-Servent, J F P J Abascal, J Aguirre, S Arridge, T Correia, J Ripoll, M Desco, J J Vaquero. J Biomed Opt, 18(7):076016, 2013.** DOI: http://dx.doi.org/10.1117/1.JBO.18.7.076016 . We use a smaller problem than the one solved in the paper to reduce computaion time to few seconds.  

##  Repository files ##

The repository contains the following files:

- **Demo_ART_SB_Reconstruction.m:** Demo to solve the 3D TV reconstruction problem for fluorescence diffuse optical tomography

- **ARTReconstruction.m:** Algebraic reconstruction technique for linear measurements

- **TV_SB_denoising_3D.m:** Solves the 3D TV denoising problem

- **TV_SB_denoising_3D.m:** Solves the 2D TV denoising problem

- **Plot2DMapsGridSolution.m:** Simple routine to display results

- **DataRed.mat:** Data, Jacobian and target solution 

If you use this code, please cite Chamorro-Servent et al. Use of Split Bregman denoising for iterative reconstruction in fluorescence diffuse optical tomography. J Biomed Opt, 18(7):076016, 2013. http://dx.doi.org/10.1117/1.JBO.18.7.076016. If you need to contact the author, please do so at juanabascal78@gmail.com, juchamser@gmail.com, desco@hggm.es

