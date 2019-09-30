
% JRM_Lite
% Copyright University College London 2015, 2016, 2017, 2018, 2019
% Author: Alexandre Bousse, Institute of Nuclear Medicine, UCL
% For research purpose only.


clear all
close all

addpath('./mex')
addpath('./mat')



N = 64 ;  NZ = 30 ; nGates = 1 ;  isTOF = 1 ;

param = jrm_makeParam(N,NZ,nGates,isTOF) ;
param.useGPU = 1 ;

[act_t,mu_t,mutilde] = jrm_makePhantom(param) ;

sino_t = jrm_makeData(act_t,mu_t,param) ;

mu = jrm_pickTimeFrame(mu_t,nGates,1) ;

f = ones(size(mu)) ;

alpha.X = zeros([param.Nspl,nGates]) ; alpha.Y = zeros([param.Nspl,nGates]) ; alpha.Z = zeros([param.Nspl,nGates]) ;


param.reinit = 1 ;
param.dispEM = 1 ;

% 
param.nEM = 30 ;
param.nIterMotion = 60 ;

% param.nEM = 2 ;
% param.nIterMotion = 2 ;


% betaQuadMotion = 1 ;
% param.betaQuadMotion.X = betaQuadMotion ; 
% param.betaQuadMotion.Y = betaQuadMotion ; 
% param.betaQuadMotion.Z = betaQuadMotion ;

param.reinitEveryNiter = 1 ;

[f,alpha] = jrm_main(f,mutilde,sino_t,alpha,param) ;
