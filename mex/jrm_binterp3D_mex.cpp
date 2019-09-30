// JRM_Lite
// Copyright University College London 2015, 2016, 2017, 2018, 2019
// Author: Alexandre Bousse, Institute of Nuclear Medicine, UCL
// For research purpose only.


#include "jrm_source_mex.h"

// mex jrm_binterp3D_mex.cpp

// warping operator

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
    // -----------------------------------------------------------------------------
    // INPUTS
    // -----------------------------------------------------------------------------
    
    
    double *mxVolForw = mxGetPr(prhs[0]) ;
    
    double *mxXnew = mxGetPr(prhs[1]) ;  // x-coordinates in the original grid
    double *mxYnew = mxGetPr(prhs[2]) ;  // y-coordinates in the original grid
    double *mxZnew = mxGetPr(prhs[3]) ;  // z-coordinates in the original grid
    
    // int *mxDimsGridBack = (int*)mxGetData(prhs[4]) ;
    //mwSize *mxDimsGridBack = (mwSize*)mxGetData(prhs[4]) ;
    mwSize *mxDimsGridBack = (mwSize*)mxGetData(prhs[4]) ;
    mwSize const *mxDimsGridForw = mxGetDimensions(prhs[1]) ; // target image dimensions (after interp)
    
    mwSize const NXback = mxDimsGridBack[1] ; 
    mwSize const NYback = mxDimsGridBack[0] ; 
    mwSize const NZback = mxDimsGridBack[2] ; 
    
    mwSize const NXforw = mxDimsGridForw[1] ; 
    mwSize const NYforw = mxDimsGridForw[0] ; 
    mwSize const NZforw = mxDimsGridForw[2] ; 
    
//     printf("NXback= %i\n",NXback) ;
//     printf("NXback= %i\n",NYback) ;
//     printf("NXback= %i\n",NZback) ;
    
    // --------------------------------------------------------------------------------
    
    // mwSize const Ngrid = mxGetM(prhs[1]) ; // size of one image edge
    

    mwSize const Ngrid_all = NXforw*NYforw*NZforw ; ;// length of mxXgrid, i.e. number of voxels
    
    int nghbCorner[3] ;
    
    int j, l, ll, lll ;
    int lz, lx, ly ;
    int nz, nx, ny ;
    double betaX, betaY, betaZ  ;
    double betaX_pc[4], betaY_pc[4], betaZ_pc[4] ; // pre-computed values
    
    double Xnew, Ynew, Znew, Xdiff, Ydiff, Zdiff ;
    double volNew, bbb ;
    
    //---------------------------------------------------
    // output: mxVolBack, i.e. W(alpha)u
    //---------------------------------------------------
    
    
    double *mxVolBack= NULL  ;
    plhs[0] = mxCreateNumericArray(3, mxDimsGridBack, mxDOUBLE_CLASS, mxREAL) ;
    mxVolBack = mxGetPr(plhs[0]) ;
    //-----------------------------------------------------
    
    for (j=0 ; j<Ngrid_all ;j++){   // can be paralellised
        
        Xnew = mxXnew[j] ; Ynew = mxYnew[j] ; Znew = mxZnew[j] ;
        initInt(nghbCorner,3) ;
        initDouble(betaX_pc,4) ; initDouble(betaY_pc,4) ; initDouble(betaZ_pc,4) ;
        findCorner(Xnew, Ynew, Znew, 1,1,1, nghbCorner) ;
        
        
        for (ll = 0; ll<4 ; ll++ ){
            nz = nghbCorner[2] + ll ;
            nx = nghbCorner[1] + ll ;
            ny = nghbCorner[0] + ll ;
            
            Zdiff = Znew - nz ;
            betaZ = beta(Zdiff) ;
            betaZ_pc[ll] = betaZ ;
            Xdiff = Xnew - nx ;
            betaX = beta(Xdiff) ;
            betaX_pc[ll] = betaX ;
            Ydiff = Ynew - ny ;
            betaY = beta(Ydiff) ;
            betaY_pc[ll] = betaY ;
            
        }
        
        
        for (lz = 0; lz<4 ; lz++){
            
            nz = nghbCorner[2] + lz ;
            if (nz>=0 && nz<NZback){
                
                betaZ = betaZ_pc[lz] ;
                for (lx = 0 ; lx<4 ; lx++  ){
                    
                    nx = nghbCorner[1] + lx ;
                    if (nx>=0 && nx<NXback){
                        
                        betaX = betaX_pc[lx] ;
                        for (ly = 0 ; ly<4 ; ly++  ){
                            ny = nghbCorner[0] + ly ;
                            
                            if (ny>=0 && ny<NYback){
                                betaY = betaY_pc[ly] ;
                                
                                lll = ny + nx*NYback +  nz*NXback*NYback ;
                                bbb = betaX * betaY * betaZ ;
                                if (bbb>0){
                                    
                                    mxVolBack[lll] = mxVolBack[lll] + mxVolForw[j]*bbb ;
                                }
                                
                            }
                            
                        }
                        
                    }
                    
                }
                
            }
            
        }
        
    }
    
}
