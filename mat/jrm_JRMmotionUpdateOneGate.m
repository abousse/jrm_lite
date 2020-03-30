% JRM_Lite
% Copyright University College London 2015, 2016, 2017, 2018, 2019
% Author: Alexandre Bousse, Institute of Nuclear Medicine, UCL
% For research purpose only.


function alpha = jrm_JRMmotionUpdateOneGate(alpha,vol,mu,sino_t,param,gateNumber)

Nspl = param.Nspl  ;
nIter = param.nIterMotion ;

alphaX = jrm_pickTimeFrame(alpha.X,param.nGates,gateNumber) ;   %
alphaY = jrm_pickTimeFrame(alpha.Y,param.nGates,gateNumber) ;
alphaZ = jrm_pickTimeFrame(alpha.Z,param.nGates,gateNumber) ;

NXspl = Nspl(1) ; NYspl = Nspl(2) ; NZspl = Nspl(3) ;
nSpline_all = NXspl*NYspl*NZspl  ;
alphaVect = zeros( 3*nSpline_all,1  ) ;

alphaVect(   1                   :     nSpline_all  ) = reshape(alphaX,nSpline_all,1) ;
alphaVect(   1  +   nSpline_all  :   2*nSpline_all  ) = reshape(alphaY,nSpline_all,1) ;
alphaVect(   1  + 2*nSpline_all  :   3*nSpline_all  ) = reshape(alphaZ,nSpline_all,1) ;


lb = -60*ones(nSpline_all*3,1) ;
ub =  60*ones(nSpline_all*3,1) ;

fopti1var = @(alphaVect)jrm_JRMmotionCostOneGate(alphaVect,vol,mu,sino_t,param,gateNumber) ;

% ------------------------------------------------------------------------------------------------------------------------
% OPTIMISATION OF THE OBJECTIVE FUNCTION WITH RESPECT TO THE MOTION PARAMETER ALPHA:

% Using 'fmincon' from the Optimization Toolbox
opt_fmincon = optimoptions(@fmincon,'Display','iter-detailed','algorithm','interior-point',...
    'GradObj','on','Hessian','lbfgs','MaxIter',nIter,'DerivativeCheck','off') ;
alphaVect = fmincon(fopti1var,alphaVect,[],[],[],[],lb,ub,[],opt_fmincon) ; % OPTIMISATION


% Alternativelly the LBFGS Toolbox by Stephen Becker can be used.
% The LBFGS toolbox MEX Wrapper should be installed and added to the path. It is available on Mathwork
% https://www.mathworks.com/matlabcentral/fileexchange/35104-lbfgsb-l-bfgs-b-mex-wrapper
% Also available on Github
% https://fr.mathworks.com/matlabcentral/fileexchange/35104-lbfgsb-l-bfgs-b-mex-wrapper
%
% opts.x0 = alphaVect ;
% opts.maxIts = nIter ;
% opts.printEvery = 1 ;
% opts.verbose = 1;
% opts.m = 9 ;
% opts.factr = 1e7 ;
% [alphaVect,~,~] = lbfgsb(fopti1var, lb, ub, opts) ; % OPTIMISATION

% ------------------------------------------------------------------------------------------------------------------------



alphaX = reshape( alphaVect(  1                   :     nSpline_all   ),  NYspl,NXspl,NZspl  ) ;
alphaY = reshape( alphaVect(  1 +   nSpline_all   :   2*nSpline_all   ),  NYspl,NXspl,NZspl  ) ;
alphaZ = reshape( alphaVect(  1 + 2*nSpline_all   :   3*nSpline_all   ),  NYspl,NXspl,NZspl  ) ;

alpha.X(:,:,:,gateNumber) = alphaX ;
alpha.Y(:,:,:,gateNumber) = alphaY ;
alpha.Z(:,:,:,gateNumber) = alphaZ ;
    




