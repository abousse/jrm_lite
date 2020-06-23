% JRM_Lite
% Copyright University College London 2015, 2016, 2017, 2018, 2019
% Author: Alexandre Bousse, Institute of Nuclear Medicine, UCL
% For research purpose only.


clear all
close all

addpath('./mex')
addpath('./mat')

 
N = 128 ;  NZ = 60 ; nGates = 3 ;  isTOF = 0 ;

param = jrm_makeParam(N,NZ,nGates,isTOF) ;
param.useGPU = 0 ; % switch to 1 if you want to use GPU
param.nEM = 40 ;
param.nIterMotion = 20 ;
param.nIterTotal = 3 ;

% activity and attenuation phantoms (3 gates)
[act_t,mu_t] = jrm_makePhantom(param) ;

% gated PET sinogram
sino_t = jrm_makeData(act_t,mu_t,param) ; 

% input mu-map for reconstruction
mu = jrm_pickTimeFrame(mu_t,nGates,1) ; 

% initialisation of the acivity image to be reconstucted
f_init = ones(size(mu)) ; 

% initialisation of the acivity image to be reconstucted
alpha.X = zeros([param.Nspl,nGates]) ; alpha.Y = zeros([param.Nspl,nGates]) ; alpha.Z = zeros([param.Nspl,nGates]) ;

% Joint image reconstruction and motion estimation
[f,alpha] = jrm_main(f_init,mu,sino_t,alpha,param) ;

% reconstruction without motion compensation (for comparison)
sino_nongated = sum(sino_t,4) ;
param.nGates = 1 ; % this is required as the gates have been summed 
f_nomc = jrm_penMLEM(f_init,sino_nongated,mu,param,1) ;

figure(3)  
imagesc(squeeze( f_nomc(:,round(size(f_nomc,2)/2),:))'  ) ;  axis image ; colormap hot ; set(gca,'XTick',[],'YTick',[]) ;
title('Reconstruction without motion compensation')






