% JRM_Lite
% Copyright University College London 2015, 2016, 2017, 2018, 2019
% Author: Alexandre Bousse, Institute of Nuclear Medicine, UCL
% For research purpose only.


function RES = pm_pickTimeFrame(A,nGates,t)

if max(t) > nGates
   error('t cannot be greater than nGates') 
end

dim = size(A) ;

if nGates == 1
    RES = A ;
    
else
    B = reshape(A,prod(dim(1:end-1)),nGates) ;
    RES = reshape( B(:,t), [dim(1:end-1),length(t)]) ;
    
end
