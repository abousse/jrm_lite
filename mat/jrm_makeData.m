% JRM_Lite
% Copyright University College London 2015, 2016, 2017, 2018, 2019
% Author: Alexandre Bousse, Institute of Nuclear Medicine, UCL
% For research purpose only.



function sino_t = jrm_makeData(act_t,mu_t,param)

N = param.dimIm(1) ;
NZ = param.dimIm(3) ;

if param.isTOF
   sino_t = zeros(floor(N/param.wTOF) + 1,N,NZ,param.nPhi,param.nGates) ;
else
   sino_t = zeros(N,NZ,param.nPhi,param.nGates) ;
end


% 
sinoSize = size(sino_t) ;
sino_t = reshape(sino_t,[] ,param.nGates) ; 

% TO DO
% chage the data generation and the rest of the code so that the background
% is proportional to the gate duration

for t = 1 : param.nGates
    attnC = jrm_makeAttnC(jrm_pickTimeFrame(mu_t,param.nGates,t), param) ;
    proj = param.gateDuration(t)*jrm_forwardProj(jrm_pickTimeFrame(act_t,param.nGates,t),param) ;
    bckg = jrm_pickTimeFrame(param.bckg,param.nGates,t) ;
    
    if param.isTOF
        attnC = reshape(attnC,[1,size(attnC)]) ;
    end
    
    expProj =    bsxfun(@times,attnC,proj)  + bckg ;
    sino_t(:,t) = expProj(:) ;
end

sino_t = reshape(sino_t, sinoSize);
