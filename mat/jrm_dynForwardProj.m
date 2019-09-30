% JRM_Lite
% Copyright University College London 2015, 2016, 2017, 2018, 2019
% Author: Alexandre Bousse, Institute of Nuclear Medicine, UCL
% For research purpose only.



function g_t = jrm_dynForwardProj(f,param,alpha)


N = size(f,1) ;
alphaX_t = alpha.X ; alphaY_t = alpha.Y ;
%nPhi = length(param.phi) ;

alphaZ_t = alpha.Z ;
NZ = size(f,3) ;

if param.isTOF
    g_t = zeros(floor(N/param.wTOF) + 1,N,NZ,param.nPhi,param.nGates) ;
else
    g_t = zeros(N,NZ,param.nPhi,param.nGates) ;
end

dimProj = size(g_t) ;
g_t = reshape(g_t,[] , param.nGates) ;

for t = 1 : param.nGates
    
    
   alphaX = alphaX_t(:,:,:,t) ; alphaY = alphaY_t(:,:,:,t) ; alphaZ = alphaZ_t(:,:,:,t) ;
   [Xnew,Ynew,Znew] = jrm_bsplineTransform_mex(param.X, param.Y, param.Z, param.Wspl,param.Nspl, alphaX, alphaY,alphaZ) ;
   f_new = jrm_fwarp3D_mex(f,Xnew,Ynew,Znew) ;
        
   proj = jrm_forwardProj(f_new,param) ;

   g_t(:,t) = proj(:) ;
    
end

g_t = reshape(g_t,dimProj) ;









