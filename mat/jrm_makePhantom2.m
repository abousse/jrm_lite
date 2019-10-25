% JRM_Lite
% Copyright University College London 2015, 2016, 2017, 2018, 2019
% Author: Alexandre Bousse, Institute of Nuclear Medicine, UCL
% For research purpose only.


function [act_t,mu_t,mutilde] = jrm_makePhantom2(param)

dimIm = param.dimIm ;
N = dimIm(1) ; NZ = dimIm(3) ;



mu_water = 0.01 ; %mm^-1

f = makeCube(N,NZ) ;
mu = makeCube(N,NZ)*mu_water*param.voxSize ; 

act_t = zeros(N,N,NZ) ;
mu_t = zeros(N,N,NZ) ;



for l = 1 : param.nGates
    
    
    shift_vect = [1,0,0]* round(  (N/5)*(l-1)/param.nGates) ;
    act_t(:,:,:,l) = circshift(f,shift_vect) ;
    mu_t(:,:,:,l) = circshift(mu,shift_vect) ;
    
end


shift_vect = [0,0,1]* round(  (NZ/8)) ;
mutilde = circshift(mu,shift_vect) ;

end





function RES = makeCube(N,NZ) 

RES = zeros(N,N,NZ) ;
RES(round(N/4):3*round(N/4) , round(N/4):3*round(N/4) ,round(NZ/3):2*round(NZ/3)  ) = 1 ;


end

