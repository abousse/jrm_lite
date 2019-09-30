% JRM_Lite
% Copyright University College London 2015, 2016, 2017, 2018, 2019
% Author: Alexandre Bousse, Institute of Nuclear Medicine, UCL
% For research purpose only.



function f_pl = jrm_mergeEMandPrior(f_em,f_reg,beta,p)

% maximum for p(f_em log(x) - x)  -   0.5*beta*( x - f_reg )^2

if (sum(beta(:))>0)
    f_pl = (   (beta.*f_reg - p) +  sqrt( (beta.*f_reg - p).^2 + 4*beta.*p.*f_em   )   )./beta/2 ;
else
    f_pl = f_em ;
end

% if p=0, then f_pl=f_reg (normally, f_em would be zero by definition)
% note that if f_em=0 and p~=0 (it happens when f_pl^old=0) and beta>0, then f_pl~=f_reg
% 

% If f_em = 0, then f_pl can be negative. This happens when the surrogate is
% of the form p*f - beta*(f-f_reg)^2 (absence of log, i.e. no counts on LOR intersecting
% the voxel). However, f_pl as calculated here (root of a 2nd order polynomial) cannot be negative

% when (beta.*f_reg - p) is negative, and 4*beta.*p.*f_em very small, then
% the numerator is similar to a + abs(a) with a<0, which returns zero
