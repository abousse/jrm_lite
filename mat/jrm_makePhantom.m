% JRM_Lite
% Copyright University College London 2015, 2016, 2017, 2018, 2019
% Author: Alexandre Bousse, Institute of Nuclear Medicine, UCL
% For research purpose only.


function [act_t,mu_t] = jrm_makePhantom(param)

dimIm = param.dimIm ;
N = dimIm(1) ; NZ = dimIm(3) ;

c1 = [round(N/2),round(N/2),round(NZ/2)] ;
c2 = [round(N/2)+round(NZ/5),round(N/2),round(NZ/2)] ;

mu_water = 0.01 ; %mm^-1

f = makeSphere(c1,NZ/3,N,NZ) ;
f = f + 4*makeSphere(c2,NZ/12,N,NZ) ;
mu = makeSphere(c1,NZ/3,N,NZ)*mu_water*param.voxSize ; 

act_t = zeros(N,N,NZ) ;
mu_t = zeros(N,N,NZ) ;



for l = 1 : param.nGates
    
    
    shift_vect = [1,0,0]* round(  (N/5)*(l-1)/param.nGates) ;
    act_t(:,:,:,l) = circshift(f,shift_vect) ;
    mu_t(:,:,:,l) = circshift(mu,shift_vect) ;
    
end



end





function RES = makeSphere(c,r,N,NZ) 

RES = zeros(N,N,NZ) ;
r2 = r*r ;
for i = 1 : N
    for j = 1 : N
        for k = 1 : NZ
            
            d2 = (i - c(1))^2 +  (j - c(2))^2 + (k - c(3))^2 ;
            
            if (d2<r2)
                RES(i,j,k) = 1 ;
            end
            
        end
    end
end


end

