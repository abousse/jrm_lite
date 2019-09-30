% JRM_Lite
% Copyright University College London 2015, 2016, 2017, 2018, 2019
% Author: Alexandre Bousse, Institute of Nuclear Medicine, UCL
% For research purpose only.

function attnC = jrm_makeAttnC(mu,param,varargin)
%attnC = pm_makeAttnC(mu,paramProj,isub,Nsub,paramRecoJRM,alpha) ;


param_mu = param ;
param_mu.isTOF = 0 ;


if nargin<3
    
    attnC = exp(-jrm_forwardProj(mu,param_mu)) ;
    
else
    
    
    alpha = varargin{1} ;
    attnC = exp(-jrm_dynForwardProj(mu,param_mu,alpha)) ;
    
end








