// JRM_Lite
// Copyright University College London 2015, 2016, 2017, 2018, 2019
// Author: Alexandre Bousse, Institute of Nuclear Medicine, UCL
// For research purpose only.

#include "jrm_source_mex.h"

// mex jrm_bsplineTransform_mex.cpp


void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
    
    double *mxXgrid = mxGetPr(prhs[0]) ;  // x-coordinates in the original grid
    double *mxYgrid = mxGetPr(prhs[1]) ;  // y-coordinates in the original grid
    double *mxZgrid = mxGetPr(prhs[2]) ;  // z-coordinates in the original grid
    
    double *mxWspl = mxGetPr(prhs[3]) ;
    int *mxNspl = (int*)mxGetData(prhs[4]) ;
    
    double *mxAlphaX = mxGetPr(prhs[5]) ;  // x bspline coefficient
    double *mxAlphaY = mxGetPr(prhs[6]) ;  // y bspline coefficient
    double *mxAlphaZ = mxGetPr(prhs[7]) ;  // z bspline coefficient
    
    mwSize const *dimsGrid = mxGetDimensions(prhs[0]) ; // image dimensions
    mwSize const NX = dimsGrid[1] ;
    mwSize const NY = dimsGrid[0] ;
    mwSize const NZ = dimsGrid[2] ;
    mwSize const Ngrid_all = NX*NY*NZ ;// length of mxXgrid, i.e. number of voxels
    
    
    mwSize const *dimsSpl = mxGetDimensions(prhs[0]) ;
    mwSize const NXspl = mxNspl[0] ;
    mwSize const NYspl = mxNspl[1] ;
    mwSize const NZspl = mxNspl[2] ;
    
    
    mwSize const Nspl_all = NXspl*NYspl*NZspl ; // size of mxXspl, i.e. number Bspline nodes
    
    int lz, lx, ly ;
    int nz, nx, ny ;
    double betaX, betaY, betaZ  ;
    double betaXZ ;
    double betaX_pc[4], betaY_pc[4], betaZ_pc[4] ; // pre-computed values
    int nghbCorner[3] ;
    
    int nghbX[64], nghbY[64], nghbZ[64] ;
    int k, l, ll, lll ;
    double X,Y,Z, Xdiff, Ydiff, Zdiff, bbb ;
    double const gridSplineWidth_X = mxWspl[0] ;
    double const gridSplineWidth_Y = mxWspl[1] ;
    double const gridSplineWidth_Z = mxWspl[2] ;
    
    double const gridSplineWidthInv_X = 1/gridSplineWidth_X ;
    double const gridSplineWidthInv_Y = 1/gridSplineWidth_Y ;
    double const gridSplineWidthInv_Z = 1/gridSplineWidth_Z ;
    
    double Xnew, Ynew, Znew ;
    
    
    //---------------------------------------------------
    // output: mxXnew, mxYnew, mxZnew
    //---------------------------------------------------
    
    mwSize const dims[3] = {NY,NX,NZ} ;
    double *mxXnew = NULL  ;
    double *mxYnew = NULL  ;
    double *mxZnew = NULL  ;
    
    plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL) ;
    plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL) ;
    plhs[2] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL) ;
    
    mxXnew = mxGetPr(plhs[0]) ;
    mxYnew = mxGetPr(plhs[1]) ;
    mxZnew = mxGetPr(plhs[2]) ;
    
    

    
    //-----------------------------------------------------
    
    for (k=0 ; k<Ngrid_all ; k++){ 
        
        X = mxXgrid[k] ; Y = mxYgrid[k] ;  Z = mxZgrid[k] ;
        initInt(nghbCorner,3) ;
        initDouble(betaX_pc,4) ; initDouble(betaY_pc,4) ; initDouble(betaZ_pc,4) ;
        findCorner(X, Y, Z, gridSplineWidth_X,gridSplineWidth_Y,gridSplineWidth_Z, nghbCorner) ;
        
        Xnew = X ; Ynew = Y ; Znew = Z ;

        // the deformation grid contains the voxel grids, plus une control point on the left and on the right of each axis. The voxel (0 0 0)
        // coincides with control point (1 1 1). The topleft neighbour of voxel (0,0,0) is the control point of coordinates 
        // (-1*gridSplineWidth_X,-1*gridSplineWidth_Y,-1*gridSplineWidth_Z) 
        
        for (ll = 0; ll<4 ; ll++ ){
            
            // Here (nx,ny,nz) starts at (-1,-1,-1).
            // In this precomputation loop we keep (nx,ny,nz) in the coordinate system centred at voxel (0,0,0), because we are only interested at their distance with (X,Y,Z). However,
            // in order to pick the correct alpha value, they will have to be shifted by +1 (see the next loop)

            nz = nghbCorner[2] + ll  ;
            nx = nghbCorner[1] + ll  ;
            ny = nghbCorner[0] + ll  ;
            
            Zdiff = Z - nz*gridSplineWidth_Z ;
            betaZ = beta(Zdiff*gridSplineWidthInv_Z) ;
            betaZ_pc[ll] = betaZ ;
            
            Xdiff = X - nx*gridSplineWidth_X ;
            betaX = beta(Xdiff*gridSplineWidthInv_X) ;
            betaX_pc[ll] = betaX ;
            
            Ydiff = Y - ny*gridSplineWidth_Y ;
            betaY = beta(Ydiff*gridSplineWidthInv_Y) ;
            betaY_pc[ll] = betaY ;
            //printf("betaY = %f\n",betaY) ;
        }

        
        
        for (lz = 0; lz<4 ; lz++){
            
            nz = nghbCorner[2] + lz + 1 ;
            if (nz>=0 && nz<NZspl){
                
                betaZ = betaZ_pc[lz] ;
                for (lx = 0 ; lx<4 ; lx++  ){
                    
                    nx = nghbCorner[1] + lx + 1 ;
                    if (nx>=0 && nx<NXspl){
                        
                        betaXZ = betaX_pc[lx]*betaZ ;
                        for (ly = 0 ; ly<4 ; ly++  ){
                            ny = nghbCorner[0] + ly + 1 ;
                            
                            if (ny>=0 && ny<NYspl){
                                bbb = betaY_pc[ly]*betaXZ ;
                                
                                lll = ny + nx*NYspl +  nz*NYspl*NXspl ;
                               
                                
                                if (bbb>0){
                                    
                                    Xnew = Xnew + mxAlphaX[lll] * bbb ;
                                    Ynew = Ynew + mxAlphaY[lll] * bbb ;
                                    Znew = Znew + mxAlphaZ[lll] * bbb ;
                                    
                                }
                                
                            }
                            
                        }
                        
                    }
                    
                }
                
            }
            
        }
        
        mxXnew[k] = Xnew ; mxYnew[k] = Ynew ; mxZnew[k] = Znew ;
        
        
    }
    
    
}









