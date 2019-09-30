% JRM_Lite
% Copyright University College London 2015, 2016, 2017, 2018, 2019
% Author: Alexandre Bousse, Institute of Nuclear Medicine, UCL
% For research purpose only.


function [f,g] = jrm_JRMmotionCostOneGate(alphaVect,vol,mu,sino_t,param,gateNumber)

% alpha corresponds to a single gate given by 'gateNumber'


param_mu = param ;
param_mu.isTOF = 0 ;

Nspl = param.Nspl  ;
Wspl = param.Wspl  ;

nGates = param.nGates ;
param.nGates = 1 ;
%blankSino = ones(size(jrm_pickTimeFrame(sino_t,nGates,1))) ;

NXspl = Nspl(1) ; NYspl = Nspl(2) ; NZspl = Nspl(3) ;
X = param.X ; Y = param.Y ; Z = param.Z ;
    
nSpline_all = NXspl*NYspl*NZspl  ;
    
alpha.X = reshape(alphaVect(   1                  :     nSpline_all  ), NYspl,NXspl,NZspl) ;
alpha.Y = reshape(alphaVect(   1  + nSpline_all   :   2*nSpline_all  ), NYspl,NXspl,NZspl) ;
alpha.Z = reshape(alphaVect(   1  + 2*nSpline_all :   3*nSpline_all  ), NYspl,NXspl,NZspl) ;
    
    
[Xnew,Ynew,Znew] = jrm_bsplineTransform_mex(X, Y, Z, Wspl, int32(Nspl), alpha.X,alpha.Y,alpha.Z) ;
    

[Www,wWw,wwW] = jrm_warp3Dderivatives_mex(vol,Xnew,Ynew,Znew) ;
Www = 1*Www ; wWw = 1*wWw ; wwW = 1*wwW ; % to remove zero values

[Www_mu,wWw_mu,wwW_mu] = jrm_warp3Dderivatives_mex(mu,Xnew,Ynew,Znew) ;
Www_mu = 1*Www_mu ; wWw_mu = 1*wWw_mu ; wwW_mu = 1*wwW_mu ; % to remove zero values


Jalpha = jrm_bsplineJac_mex(X,Y,Z,Wspl, Nspl) ;
dimSino_t = size(sino_t) ;

if nGates == 1  % can happen often!
    dimSino_t = [dimSino_t,1] ;
end
% doesn't deal with single angle though....


if param.isTOF
    dimSinoNoTOF = [1,dimSino_t(2:end-1)] ;
else
    dimSinoNoTOF = dimSino_t(1:end-1) ;
end


bckg = jrm_pickTimeFrame(param.bckg,nGates,gateNumber) ;
sino = jrm_pickTimeFrame(sino_t,nGates,gateNumber) ;
attnC = jrm_makeAttnC(mu,param,alpha) ;
attnC = reshape(attnC,dimSinoNoTOF) ;


ACMC_proj = param.gateDuration(gateNumber)*bsxfun(@times,...
    jrm_dynForwardProj(vol,param,alpha),attnC) ;

f_likelihood = sino.*log(ACMC_proj + bckg ) - ACMC_proj ;
f_likelihood(isnan(f_likelihood)) = 0 ; % because 0*log(0) = NaN
f_likelihood = - sum(f_likelihood(:)) ; % negative likelihood!

V = sino./ (ACMC_proj + bckg) - 1 ;

V(isnan(V)) = -1 ;

% -----------------------------------------------------------------------------------------------
grad_f = param.gateDuration(gateNumber)*jrm_backProj(...
    bsxfun(@times,V,attnC),param) ;
grad_f = grad_f(:) ;



if param.isTOF
    grad_mu = jrm_backProj(squeeze(sum(ACMC_proj.*V,1)),param_mu) ; grad_mu = grad_mu(:) ;
else
    grad_mu = jrm_backProj(ACMC_proj.*V,param_mu) ; grad_mu = grad_mu(:) ;
end

clear ACMC_proj
clear V

gradMotion_X  = Www*grad_f  - Www_mu*grad_mu ;
g_likelihood_X = Jalpha*gradMotion_X ; clear gradMotion_X
gradMotion_Y  = wWw*grad_f  - wWw_mu*grad_mu ;
g_likelihood_Y = Jalpha*gradMotion_Y ; clear gradMotion_Y
gradMotion_Z  = wwW*grad_f  - wwW_mu*grad_mu ;
g_likelihood_Z = Jalpha*gradMotion_Z ; clear gradMotion_Z


% negative likelihood
g_likelihood(1                :  nSpline_all  ) = -g_likelihood_X ;  clear g_likelihood_X
g_likelihood(nSpline_all+1    :  2*nSpline_all) = -g_likelihood_Y ;  clear g_likelihood_Y
g_likelihood(2*nSpline_all+1  :  3*nSpline_all) = -g_likelihood_Z ;  clear g_likelihood_Z


clear grad_mu
clear grad_f

f =  f_likelihood ;
g =  g_likelihood ;



% ------------------------------------------------------------------------------------------------------
% Quadratic prior
% ------------------------------------------------------------------------------------------------------

[f_quad_X , g_quad_X] = jrm_quadPrior(alpha.X) ;
[f_quad_Y , g_quad_Y] = jrm_quadPrior(alpha.Y) ;
[f_quad_Z , g_quad_Z] = jrm_quadPrior(alpha.Z) ;

f_quad =   param.betaQuadMotion.X*f_quad_X + param.betaQuadMotion.Y*f_quad_Y + param.betaQuadMotion.Z*f_quad_Z   ;
g_quad = [ param.betaQuadMotion.X*g_quad_X ; param.betaQuadMotion.Y*g_quad_Y ; param.betaQuadMotion.Z*g_quad_Z ] ;



f = f + param.gateDuration(gateNumber)*f_quad ;
g = g(:) + param.gateDuration(gateNumber)*g_quad(:) ;















