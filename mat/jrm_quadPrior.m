% JRM_Lite
% Copyright University College London 2015, 2016, 2017, 2018, 2019
% Author: Alexandre Bousse, Institute of Nuclear Medicine, UCL
% For research purpose only.




function [f,g, im_reg,omega_sum] = jrm_quadPrior(im,varargin)

% warning :
% here f is the cost function value and g is the gradient

% compute quadratic smoothed version of a volume, in order to be used for
% separable surrogate



list = [0 1 -1] ;



vect = zeros(27,3) ;
ii = 0 ;

for i = 1 : 3
    for j = 1 : 3
        for k = 1 : 3
            ii = ii + 1 ;
            vect(ii,:) = [list(i) list(j) list(k)] ;
        end
    end
end

vect(1,:) = [] ;
omega_sum = 0 ;

F = 0 ;
im_av = 0 ;
im_reg = 0 ;

for i = 1 : size(vect,1)
    
    shift_vect = vect(i,:) ;
    im_shifted = circshift(im,shift_vect) ;
    
    omega = ones(size(im))/norm(shift_vect) ; 
    
   
    F = F + 0.5 * omega .* (im - im_shifted).^2 ;
    omega_sum = omega_sum + omega ;
    
    im_av = im_av + im_shifted.*omega ;
    im_reg = im_reg + (im + im_shifted).*omega ;
    
end

im_av = im_av./omega_sum ;
im_reg = 0.5*im_reg./omega_sum ;

f = sum(F(:)) ;
g = 2*omega_sum .* (im - im_av) ; % the 2 is correct!
g = g(:) ;










