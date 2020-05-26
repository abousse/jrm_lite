% JRM_Lite
% Copyright University College London 2015, 2016, 2017, 2018, 2019
% Author: Alexandre Bousse, Institute of Nuclear Medicine, UCL
% For research purpose only.



function im = jrm_backProj(sino,param)
%im = pm_backProj(sino,param,isub,Nsub)
% General back projector

% 
if param.isTOF
    N = size(sino,2) ;
else
    N = size(sino,1) ;
end

sigma_mm = param.FWHM/2.3555 ;
sigma_time_mm = param.FWHM_TOF/2.3555 ;
sigma_vox = sigma_mm/param.voxSize ;
sigma_time_vox = sigma_time_mm/param.voxSize ;
phi = 180*param.phi/pi ;
nPhi = length(phi) ; 

[I,J] = ndgrid(1:N,1:N);
ic = floor((N+1)/2)  ;
jc = ic ;
D = double((I-ic).^2+(J-jc).^2<=round(N/2)^2);

nPSF = round(N/5) ;
nPSF_time = N ;

if (mod(nPSF,2)==0)
    h = fspecial('gaussian',[1 nPSF-1],sigma_vox) ;
else
    h = fspecial('gaussian',[1 nPSF],sigma_vox) ;
end

if (mod(nPSF_time,2)==0)
    h_time = fspecial('gaussian',[1 nPSF_time-1],sigma_time_vox) ;
else
    h_time = fspecial('gaussian',[1 nPSF_time],sigma_time_vox) ;
end



if param.isTOF
    NZ = size(sino,3) ;
else
    NZ = size(sino,2) ;
end

H = {h,h,h} ;
H_time = {h_time,1,1} ;

im = zeros(N,N,NZ) ;

if (param.useGPU == 1)
    im = gpuArray(single(im)) ;
    if param.isTOF~=1
        sino = gpuArray(single(sino)) ;
    end
end

for i = 1 : nPhi
    
    if param.isTOF
        r_sino = jrm_oneDimBbinning_mex(sino(:,:,:,i),param.wTOF,N) ;%/param.wTOF ;
        
        if (param.useGPU == 1)
            r_sino = gpuArray(single(r_sino))  ;
        end
        r_sino = convnsep(H_time,r_sino,'same') ;
    else
        r_sino = repmat(reshape(sino(:,:,i),1,N,NZ  )   ,[N 1 1]) ;
    end
    
    
    im = im + imrotate(r_sino,-phi(i),'bilinear','crop') ;
    
end

im = im.*D ;
im = convnsep(H,im,'same') ;

if (param.useGPU == 1)
    im = double(gather(im)) ;
else
    im = double(im) ;
end


    
    
    
    
    
