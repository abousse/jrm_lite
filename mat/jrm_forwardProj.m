% JRM_Lite
% Copyright University College London 2015, 2016, 2017, 2018, 2019
% Author: Alexandre Bousse, Institute of Nuclear Medicine, UCL
% For research purpose only.




function sino = jrm_forwardProj(im,paramProj)
% General forward projector






N = size(im,1) ;
sigma_mm = paramProj.FWHM/2.3555 ;
sigma_time_mm = paramProj.FWHM_TOF/2.3555 ;
sigma_vox = sigma_mm/paramProj.voxSize ;
sigma_time_vox = sigma_time_mm/paramProj.voxSize ;
phi = 180*paramProj.phi/pi ;

nPhi = length(phi) ; % now nPhi is the size of the subset....

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



NZ = size(im,3) ;
H = {h,h,h} ;
H_time = {h_time,1,1} ;

if paramProj.isTOF
    %sino = zeros(N,N,NZ,nPhi) ;
    sino = zeros(floor(N/paramProj.wTOF) + 1, N,NZ,nPhi) ; 
else
    sino = zeros(N,NZ,nPhi) ;
end

if (paramProj.useGPU == 1)
    im = gpuArray(single(im)) ;
    sino = gpuArray(single(sino)) ;
end

im = convnsep(H,im,'same') ;


for i = 1 : nPhi
    
    im_rot = imrotate(im,phi(i),'bilinear','crop') ;
    
    if paramProj.isTOF
        im_rot_conv = double(gather(convnsep(H_time,im_rot,'same'))) ;
        sino(:,:,:,i) = jrm_oneDimFbinning_mex(im_rot_conv,paramProj.wTOF) ;%/paramProj.wTOF ;
    else
        sino(:,:,i) = squeeze(sum(im_rot,1)) ;
    end
    
end

if paramProj.useGPU
    sino = double(gather(sino)) ;
end

    




    
