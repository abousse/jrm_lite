% JRM_Lite
% Copyright University College London 2015, 2016, 2017, 2018, 2019
% Author: Alexandre Bousse, Institute of Nuclear Medicine, UCL
% For research purpose only.



function param = jrm_makeParam(N,NZ,nGates,isTOF)


%% Sinogram and projection parameters 

param.isTOF = isTOF ;
param.dimIm = [N,N,NZ] ;
param.voxSize = 3.125*128/N ;
param.nGates = nGates ;
param.useGPU = 0 ;
param.FWHM = 4 ; %mm
param.FWHM_TOF = 75 ; %mm
%param.FWHM_TOF = 40 ; %mm
param.phi = 0:0.04:pi ;
param.nPhi = length(param.phi) ;
param.voxSize = 3.125 ; %mm

%param.normSino = ones([param.dimIm(1),param.dimIm(3),param.nPhi]) ;


param.wTOF = 4 ; % Time of flight bin width (in voxels)
param.nt =   floor(N/param.wTOF) + 1 ;

if param.isTOF
    param.bckg = ones([param.nt,param.dimIm(1),param.dimIm(3),param.nPhi,param.nGates])*0.1 / param.nt ;
else
    param.bckg = ones([param.dimIm(1),param.dimIm(3),param.nPhi,param.nGates])*0.1  ;
end


param.gateDuration = ones(1,nGates) ;

%% Motion grid and Bsplines
width = 6 ;
[X, Y, Z] = meshgrid(0:N-1,0:N-1,0:NZ-1) ;
param.X = X ;
param.Y = Y ;
param.Z = Z ;

WsplX = width ; WsplY = WsplX ; WsplZ = WsplX ;
    
% previous values prior to the addition of additional control points
%NsplX = floor(N/WsplX) + 1 ;  NsplY = floor(N/WsplY) + 1 ;  NsplZ = floor(NZ/WsplZ) + 1 ;
    
NsplX = floor(N/WsplX) ;  NsplY = floor(N/WsplY)  ;  NsplZ = floor(NZ/WsplZ)  ;
% add one on the left, 2 on the right 
NsplX = NsplX + 3 ;
NsplY = NsplY + 3 ;
NsplZ = NsplZ + 3 ;
    
    
Wspl = [WsplX WsplY WsplZ] ; Nspl = [int32(NsplX) int32(NsplY) int32(NsplZ)] ;
param.Nspl = Nspl ;
param.Wspl = Wspl ;


%% Reconstruction parameters

% motion quadratic spatial prior weigth
betaQuadMotion = 0.01 ;
param.betaQuadMotion.X = betaQuadMotion ; 
param.betaQuadMotion.Y = betaQuadMotion ; 
param.betaQuadMotion.Z = betaQuadMotion ;

% image quadratic spatial prior weigth
param.betaQuadImage = 0.00001 ;


param.nIterMotion = 50 ; % sino-based registration inner loop

param.nEM = 40 ; 
param.nIterTotal = 10 ; % total iter
param.reinitEveryNiter = 5 ;

param.reinit = 0 ;
param.reinitFinalIter = 1 ;
param.gateNumber4init = 1 ;









