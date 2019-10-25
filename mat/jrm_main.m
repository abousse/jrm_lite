% JRM_Lite
% Copyright University College London 2015, 2016, 2017, 2018, 2019
% Author: Alexandre Bousse, Institute of Nuclear Medicine, UCL
% For research purpose only.


function [f,alpha] = jrm_main(f,mu,sino_t,alpha,param)

X = param.X ;
Y = param.Y ;
Z = param.Z ;

nGates = param.nGates ;
Nspl = param.Nspl  ;
Wspl = param.Wspl ;

%NZ = size(f,3) ;

if param.reinit == 1 
    f = ones(size(f)) ;
    disp('initialisation...')
    f = jrm_penMLEM(f,sino_t,mu,param.gateNumber4init,param,param.nEM) ;
end



for n = 1 : param.nIterTotal
    
    % =======================================================================
    % maximisation w.r.t. alpha
    % =======================================================================
    
    disp(['iteration ',num2str(n),'/',num2str(param.nIterTotal)])
    tic
    
    disp('motion estimation...')
    
    
    nIterMotion = param.nIterMotion ;
    
    for t = 1 : param.nGates
        
        disp(['gate ',num2str(t),'/',num2str(param.nGates)])
        
        if ( n==1 && t>1 )
            
            alpha.X(:,:,:,t) = alpha.X(:,:,:,t-1) ;
            alpha.Y(:,:,:,t) = alpha.Y(:,:,:,t-1) ;
            alpha.Z(:,:,:,t) = alpha.Z(:,:,:,t-1) ;
        end
        
        alpha = jrm_JRMmotionUpdateOneGate(alpha,f,mu,sino_t,param,t,nIterMotion) ;
        
    end
    
    
    % =======================================================================
    % motion-compensated reconstruction
    % =======================================================================
    
    alphaX_t = alpha.X ;
    alphaY_t = alpha.Y ;
    alphaZ_t = alpha.Z ;
    
    
    
    
    
    % test if f should be reinitialised
    if (mod(n,param.reinitEveryNiter) == 0 ||  (mod(n,param.nIterTotal) == 0  && param.reinitFinalIter == 1) || n==1)
        
        f = ones(size(f)) ;
    end
    
    disp('Motion compensated reco...')
   
    
    
    f = jrm_dynPenMLEM(f,sino_t,mu,param,alpha,param.nEM) ;
    
   % warp the new f
   
   
   for t = 1 : nGates
       
       [Xnew,Ynew,Znew] = jrm_bsplineTransform_mex(X, Y, Z, Wspl,Nspl, alphaX_t(:,:,:,t),alphaY_t(:,:,:,t),alphaZ_t(:,:,:,t)) ;
       
       Wf = jrm_fwarp3D_mex(f,Xnew,Ynew,Znew) ;
       Wmu = jrm_fwarp3D_mex(mu,Xnew,Ynew,Znew) ;
       
       figure(2)
       im1 = subplot(nGates,2,(t-1)*2 + 1) ;
       imagesc(   squeeze(    Wf(:,round(size(f,2)/2),:))'  ) ; axis(im1,'image') ; colormap(im1,'hot') ; set(gca,'XTick',[],'YTick',[]) ;
       title(['Wf, gate ',num2str(t)])
       
       
       im2 = subplot(nGates,2,(t-1)*2 + 2) ; 
       imagesc(   squeeze(    Wmu(:,round(size(f,2)/2),:))'    ) ; axis(im2,'image') ; colormap(im2,1-gray) ; set(gca,'XTick',[],'YTick',[]) ;
       title(['Wmu, gate ',num2str(t)])
       
       pause(0.1)
       
   end
   
   
   
    
end










