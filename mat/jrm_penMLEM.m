% JRM_Lite
% Copyright University College London 2015, 2016, 2017, 2018, 2019
% Author: Alexandre Bousse, Institute of Nuclear Medicine, UCL
% For research purpose only.


function f = jrm_penMLEM(f,sino_t,mu,param,gateNumber)

sino = jrm_pickTimeFrame(sino_t,param.nGates,gateNumber) ;
dimSino = size(sino) ;
%normSino = param.normSino ;

if param.nPhi == 1  % just it case but it shouldn't really happen
    dimSino = [dimSino,1] ;
end

if param.isTOF
    dimSinoNoTOF = [1,dimSino(2:end)] ;
else
    dimSinoNoTOF = dimSino;
end

attnC = jrm_makeAttnC(mu,param) ;
attnC = reshape(attnC,dimSinoNoTOF) ; % add an extra dimension of size 1 if isTOF=1

bckg = jrm_pickTimeFrame(param.bckg,param.nGates,gateNumber) ;

if param.isTOF
    
    norm_EM = jrm_backProj(...
        param.gateDuration(gateNumber)*repmat(attnC,...
        [param.nt, ones(1,length(dimSino)-1)]) , ...
        param) ;
    
else
    norm_EM = jrm_backProj(param.gateDuration(gateNumber)*attnC,param) ;
end

for k = 1 : param.nEM
    
    
    disp(['Modified EM, iter ',num2str(k),'/',num2str(param.nEM)] ) ;
    
    % expected projection
    A = jrm_forwardProj(f,param) ;
    A = param.gateDuration(gateNumber)*bsxfun(@times,A,attnC) ;
    
    A = A + bckg ;
    
    % ratio and backprojection
    A = sino./A ;
    A(isnan(A)) = 0 ;
    A(isinf(A)) = 0 ;
    
    A = param.gateDuration(gateNumber)*bsxfun(@times,A,attnC) ;
    A = jrm_backProj(A,param) ;
    
    % EM update
    A = A./norm_EM ;
    A(isnan(A)) = 0 ;
    A(isinf(A)) = 0 ;
    f_em = f.*A ; clear A
    
    
    %f_old = f ;
    
    if (param.betaQuadImage>0 && k>1)
        
        
        [~,~,f_reg,omega_sum] = jrm_quadPrior(f) ;
        f = jrm_mergeEMandPrior(f_em,f_reg,4*param.gateDuration(gateNumber)*omega_sum*param.betaQuadImage,norm_EM) ;
        
        % Here we proceeded as in Wang and Qi, IEEE TMI 2012
        % beta should be multiplied by 4 because the surrogate is
        % 2*beta*(f - f_reg)^2 and jrm_mergeEMandPrior uses 0.5*beta*(f - f_reg)^2
        % we also multiplied by param.gateDuration(gateNumber)to smooth
        % proportionally to the total acquisition time
   
    else
        f = f_em ;
    end
    
    

    
    
end



