JRM_Lite 
Copyright University College London 2015, 2016, 2017, 2018, 2019 
Author: Alexandre
Bousse, Institute of Nuclear Medicine, UCL 
Email: alexandre@bousse.fr 
For research purpose only.

This is JRM_lite, a light version of Joint Reconstruction/Motion (JRM) toolbox developed for [1], [2] 
and to some extend [3]. This code can be used for basic simulations.

There is no userguide for now so users should utilise the 2 example scripts and dig in the source code.

The XCAT phantom was replaced with a basic spherical (Death Star like) phantom due to copyrights.

FUNCTIONS THAT ARE NOT INCLUDED IN THIS PACKAGE:
(i) REQUIRED: the 'convnsep.m' function by Igor Solovey for separable kernel convolution, available on Mathwork at: https://www.mathworks.com/matlabcentral/fileexchange/27957-separable-n-dimensional-convolution, must be in the path. 
(ii) OPTIONAL: the LBFGS toolbox MEX Wrapper by Stephen Becker, available on Mathwork: https://www.mathworks.com/matlabcentral/fileexchange/35104-lbfgsb-l-bfgs-b-mex-wrapper;  also available on Github: https://github.com/stephenbeckr/L-BFGS-B-C. The L-BFGS-B algorithm was written by Ciyou Zhu, R. H. Byrd, P. Lu and J. Nocedal [4]; please cite their work if you use it. The LBFGS function is used in 'jrm_JRMmotionUpdateOneGate.m' for motion estimation, but can also be replaced by 'fmincon' (MATLAB Optimization Toolbox). 

The MEX can be compiled with the compile.m script in the ./mex folder.

This code can be used by anyone for research purpose, provided the work from [1,2] is acknowledged.

scriptTest.m is a basic example of motion estimation/compensation from gated PET/CT. The data consist of a simple phantom in translation. Each of the 3 gates correspond to one position of the phantom.  

JRM_lite can utilise GPU acceleration, mostly for the rotation-based PET/CT projector (CUDA and Matlab Parallel Computing Toolbox required), and can be disabled by setting param.useGPU=1 in the script. 



References

[1] A. Bousse, O. Bertolli, D. Atkinson, S. Arridge, S. Ourselin, B. F. Hutton, and K. Thielemans. “Maximum-Likelihood Joint Image Reconstruction/Motion Estimation in Attenuation-Corrected Respiratory Gated PET/CT using a Single Attenuation Map”. In: IEEE Trans. Med. Imag. 35.1 (2016), pp. 217–228.

[2] A. Bousse, O. Bertolli, D. Atkinson, S. Arridge, S. Ourselin, B. F. Hutton, and K. Thielemans. “Maximum-Likelihood Joint Image Reconstruction and Motion Estimation with Misaligned Attenuation in TOF-PET/CT”. 
In: Phys. Med. Biol. 61.3 (2016), pp. L11–19.

[3] A. Bousse, R. Manber, B. F. Holman, D. Atkinson, S. Arridge, S. Ourselin, B. F. Hutton, and K. Thielemans. “Evaluation of a Direct Motion Estimation/Correction Method in Respiratory-Gated PET/MRI with Motion-Adjusted Attenuation”. In: Med. Phys. 44.6 (2017), pp. 2379–2390.

[4] C. Zhu, R. H. Byrd, P. Lu, and J. Nocedal. “Algorithm 778: L-BFGS-B: Fortran subroutines for large-scale bound-constrained optimization”. In: ACM T. Math. Software, 23.4 (1997), pp. 550–560.
