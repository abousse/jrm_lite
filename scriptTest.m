
% JRM_Lite
% Copyright University College London 2015, 2016, 2017, 2018, 2019
% Author: Alexandre Bousse, Institute of Nuclear Medicine, UCL
% For research purpose only.


clear all
close all

addpath('./mex')
addpath('./mat')



N = 64 ;  NZ = 30 ; nGates = 3 ;  isTOF = 0 ;

param = jrm_makeParam(N,NZ,nGates,isTOF) ;
param.useGPU = 1 ;

[act_t,mu_t] = jrm_makePhantom(param) ;

sino_t = jrm_makeData(act_t,mu_t,param) ;

mu = jrm_pickTimeFrame(mu_t,nGates,1) ;

f = ones(size(mu)) ;

alpha.X = zeros([param.Nspl,nGates]) ; alpha.Y = zeros([param.Nspl,nGates]) ; alpha.Z = zeros([param.Nspl,nGates]) ;


param.reinit = 1 ;
param.dispEM = 1 ;


param.nEM = 40 ;
param.nIterMotion = 20 ;

[f,alpha] = jrm_main(f,mu,sino_t,alpha,param) ;




