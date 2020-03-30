% JRM_Lite
% Copyright University College London 2015, 2016, 2017, 2018, 2019
% Author: Alexandre Bousse, Institute of Nuclear Medicine, UCL
% For research purpose only.


function f = jrm_dynPenMLEM(f,sino_t,mu,alpha,param)


dimSino_t = size(sino_t) ;


if param.nGates == 1  % can happen often!
    dimSino_t = [dimSino_t,1] ;
end
% doesn't deal with single angle though....


if param.isTOF
    dimSinoNoTOF = [1,dimSino_t(2:end)] ; %tof dimension removed 
else
    dimSinoNoTOF = dimSino_t ;
end


gateDuration = reshape(param.gateDuration,[ones(1,length(dimSino_t(1:end-1))),param.nGates]   ) ;

attnC_t = jrm_makeAttnC(mu,param,alpha) ;
attnC_t = reshape(attnC_t,dimSinoNoTOF) ;
bckg_t = param.bckg ;
    
if param.isTOF
    
    norm_t_EM = jrm_dynBackProj(...
        bsxfun(@times,...
        repmat(attnC_t,[param.nt, ones(1,length(dimSino_t)-1)]),...
        gateDuration),...
        param,alpha) ;
    
 else
    norm_t_EM = jrm_dynBackProj(...
        bsxfun(@times,attnC_t,gateDuration),...
        param,alpha) ;
    
 end

    
for k = 1 : param.nEM
     
     disp(['Modified Motion compensated EM, iter ',num2str(k),'/',num2str(param.nEM)]) ;
     
    % expected projection
     A = jrm_dynForwardProj(f,param,alpha) ;
     A = bsxfun(@times, A, attnC_t) ;
     A = bsxfun(@times, A, gateDuration) ;
     A = A + bckg_t ;
     
     % ratio and backprojection
     A = sino_t./A ;
     A(isnan(A)) = 0 ;
     A(isinf(A)) = 0 ;
     
     A = bsxfun(@times,A,attnC_t) ;
     A = bsxfun(@times,A,gateDuration) ;
     A = jrm_dynBackProj(A,param,alpha) ;
     
     A = A./norm_t_EM ;
     A(isnan(A)) = 0 ;
     A(isinf(A)) = 0 ;
     f_em = f.*A ; clear A
     
     if (param.betaQuadImage>0 && k>1)
         
          [~,~,f_reg,omega_sum] = jrm_quadPrior(f) ;
         f = jrm_mergeEMandPrior(f_em,f_reg,4*sum(param.gateDuration(:))*omega_sum*param.betaQuadImage,norm_t_EM) ;
         % Here we proceeded as in Wang and Qi, IEEE TMI 2012.
         % beta should be multiplied by 4 because the surrogate is 2*beta*(f - f_reg)^2 and jrm_mergeEMandPrior uses 0.5*beta*(f - f_reg)^2
         % we also multiplied by sum(param.fracGates(:)) to smooth proportionally to the total acquisition time
     else
         f = f_em ;
     end
     

     
     
     
     
 end




