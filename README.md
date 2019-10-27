JRM_Lite 
Copyright University College London 2015, 2016, 2017, 2018, 2019 
Author: Alexandre
Bousse, Institute of Nuclear Medicine, UCL 
Email: alexandre@bousse.fr 
For research purpose only.

This is JRM_lite, a light version of Joint Reconstruction/Motion (JRM) toolbox developed for [1], [2] 
and to some extend [3]. This code can be used for basic simulations.

There is no useguide for now so users should utilise the 2 example scripts and dig in the source code.

scriptTest.m is a basic example of motion estimation/compensation from gated PET/CT, where
the phantom is translated across the 3 gates, and corresponds to [1].

scriptTest2.m uses an attenuation map mutilde that differs from that of the simulated data, in order
to test activity/attenuation realignment in TOF-PET/CT. It is similar to the experiment proposed
in [2].

The XCAT phantom was replaced with a basic spherical phantom due to copyrights.
There are 2 softwares that should be downloaded: 
(i) the convnsep function for separable kernel convolution, available on Mathwork at:
https://www.mathworks.com/matlabcentral/fileexchange/27957-separable-n-dimensional-convolution 
(ii) the LBFGS toolbox MEX Wrapper, also available on Mathwork
https://www.mathworks.com/matlabcentral/fileexchange/35104-lbfgsb-l-bfgs-b-mex-wrapper.

The LBFGS function is used in jrm_JRMmotionUpdateOneGate for motion estimation, but can
also be replaced by fmincon (Optimization Toolbox).

Please do not forget to compile the MEX files with the compile.m script in the ./mex folder.

This code can be used by anyone for research purpose, provided the work from [1,2] is acknowledged.



References

[1] A. Bousse, O. Bertolli, D. Atkinson, S. Arridge, S. Ourselin, B. F. Hutton, and K. Thielemans.
“Maximum-Likelihood Joint Image Reconstruction/Motion Estimation in Attenuation-Corrected
Respiratory Gated PET/CT using a Single Attenuation Map”. 
In: IEEE Trans. Med. Imag. 35.1 (2016), pp. 217–228.

[2] A. Bousse, O. Bertolli, D. Atkinson, S. Arridge, S. Ourselin, B. F. Hutton, and K. Thielemans.
“Maximum-Likelihood Joint Image Reconstruction and Motion Estimation with Misaligned Attenuation in TOF-PET/CT”. 
In: Phys. Med. Biol. 61.3 (2016), pp. L11–19.

[3] A. Bousse, R. Manber, B. F. Holman, D. Atkinson, S. Arridge, S. Ourselin, B. F. Hutton, and K.
Thielemans. “Evaluation of a Direct Motion Estimation/Correction Method in Respiratory-Gated PET/MRI with Motion-Adjusted Attenuation”. 
In: Med. Phys. 44.6 (2017), pp. 2379–2390.
