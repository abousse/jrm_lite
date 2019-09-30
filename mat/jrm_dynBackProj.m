% JRM_Lite
% Copyright University College London 2015, 2016, 2017, 2018, 2019
% Author: Alexandre Bousse, Institute of Nuclear Medicine, UCL
% For research purpose only.


function f = jrm_dynBackProj(g_t,param,alpha)
        



alphaX_t = alpha.X ; alphaY_t = alpha.Y ; alphaZ_t = alpha.Z ;


f = 0 ;

for t = 1 : param.nGates
    
    f_t = jrm_backProj(jrm_pickTimeFrame(g_t,param.nGates,t),param) ;
    
    alphaX = alphaX_t(:,:,:,t) ; alphaY = alphaY_t(:,:,:,t) ; alphaZ = alphaZ_t(:,:,:,t) ;
    [Xnew,Ynew,Znew] = jrm_bsplineTransform_mex(param.X, param.Y, param.Z, param.Wspl,param.Nspl, alphaX, alphaY,alphaZ) ;
    f = f + jrm_bwarp3D_mex(f_t,Xnew,Ynew,Znew) ;
        
end
