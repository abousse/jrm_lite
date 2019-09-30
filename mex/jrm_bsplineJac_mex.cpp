// JRM_Lite
// Copyright University College London 2015, 2016, 2017, 2018, 2019
// Author: Alexandre Bousse, Institute of Nuclear Medicine, UCL
// For research purpose only.

#include "jrm_source_mex.h"

// mex  -largeArrayDims jrm_bsplineJac_mex.cpp

//  Bspline deformation model Jacobian (transposed)
//  it is independant of alpha (the model is linear in alpha)

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
    // -----------------------------------------------------------------------------
    // INPUTS
    // -----------------------------------------------------------------------------
    
    double *mxXgrid = mxGetPr(prhs[0]) ;  // x-coordinates in the original grid
    double *mxYgrid = mxGetPr(prhs[1]) ;  // y-coordinates in the original grid
    double *mxZgrid = mxGetPr(prhs[2]) ;  // z-coordinates in the original grid
    
    double *mxWspl = mxGetPr(prhs[3]) ;
    int *mxNspl = (int*)mxGetData(prhs[4]) ;
    
    mwSize const *dimsGrid = mxGetDimensions(prhs[0]) ; // image dimensions
    mwSize const NX = dimsGrid[1] ;
    mwSize const NY = dimsGrid[0] ;
    mwSize const NZ = dimsGrid[2] ;
    mwSize const Ngrid_all = NX*NY*NZ ;// length of mxXgrid, i.e. number of voxels
    
    
    mwSize const NXspl = mxNspl[0] ;    
    mwSize const NYspl = mxNspl[1] ;
    mwSize const NZspl = mxNspl[2] ;
    
    
    mwSize const Nspl_all = NXspl*NYspl*NZspl ; // size of mxXspl, i.e. number Bspline nodes
    
    int nghbCorner[3] ; // this variable contains the coordinates the current voxel neighbourhood corner point,
    // and is used as a starting point to find all the neighbours. In this function the neighbours are control points (not voxels)
    // and the coordinates are in [0,1,...,Nspl] x [0,1,...,Nspl] x [0,1,...,Nspl].
    int  k, ll, lll ;
    double X,Y,Z, Xdiff, Ydiff, Zdiff, bbb ;
    
    double const gridSplineWidth_X = mxWspl[0] ;
    double const gridSplineWidth_Y = mxWspl[1] ;
    double const gridSplineWidth_Z = mxWspl[2] ;
    
    double const gridSplineWidthInv_X = 1/gridSplineWidth_X ;
    double const gridSplineWidthInv_Y = 1/gridSplineWidth_Y ;
    double const gridSplineWidthInv_Z = 1/gridSplineWidth_Z ;
//    printf("gridSplineWidth= %f\n",gridSplineWidth) ;
    
    int lz, lx, ly ;
    int nz, nx, ny ;
    double betaX, betaY, betaZ  ;
    double betaX_pc[4], betaY_pc[4], betaZ_pc[4] ; // pre-computed values
    
    
    // -----------------------------------------------------------------------------
    // OUTPUT
    // -----------------------------------------------------------------------------
    double *mxSplineGrad = NULL  ;
    
    mwIndex *ir = NULL ;
    mwIndex *jc = NULL ;
    mwIndex entry_counter = 0 ; // current matrix entry address.
    
    mwSize const Nelts = Ngrid_all * 64 ;
    plhs[0] = mxCreateSparse(Nspl_all,Ngrid_all,Nelts, mxREAL) ; // Pointer to the output sparse matrix structre
    mxSplineGrad = mxGetPr(plhs[0]) ; // pointer to the first matrix element
    ir = mxGetIr(plhs[0]) ; // row indices
    jc = mxGetJc(plhs[0]) ; // colomun indices. 
    // Please read Matlab documentation for sparse matrices indexation because it is not intuitive at all
    
    // -----------------------------------------------------------------------------
    
    for (k=0 ; k<Ngrid_all ; k++){  // CANNOT be parallelised

        X = mxXgrid[k] ; Y = mxYgrid[k] ;  Z = mxZgrid[k] ; 
        
        initInt(nghbCorner,3) ;
        initDouble(betaX_pc,4) ; initDouble(betaY_pc,4) ; initDouble(betaZ_pc,4) ;
        findCorner(X, Y, Z, gridSplineWidth_X,gridSplineWidth_Y,gridSplineWidth_Z, nghbCorner) ;
        
        // pre-computation of the products beta((x_k-xtilde_l) / w) * beta((y_k-ytilde_l) / w) * beta((z_k-ztilde_l) / w)
        for (ll = 0; ll<4 ; ll++ ){
            nz = nghbCorner[2] + ll ;
            nx = nghbCorner[1] + ll ;
            ny = nghbCorner[0] + ll ;
            
            Zdiff = Z - nz*gridSplineWidth_Z ;
            betaZ = beta(Zdiff*gridSplineWidthInv_Z) ;
            betaZ_pc[ll] = betaZ ;
            
            Xdiff = X - nx*gridSplineWidth_X ;
            betaX = beta(Xdiff*gridSplineWidthInv_X) ;
            betaX_pc[ll] = betaX ;
            
            Ydiff = Y - ny*gridSplineWidth_Y ;
            betaY = beta(Ydiff*gridSplineWidthInv_Y) ;
            betaY_pc[ll] = betaY ;
        }
        
        ///////////////////////////////////////////////////////////////////
        jc[(mwIndex)(k)] = entry_counter ;
        ///////////////////////////////////////////////////////////////////     
        
        // loop in the neighbourhood, starting from corner
        // note that the axis order is y, x, z
        for (lz = 0; lz<4 ; lz++){
            
            nz = nghbCorner[2] + lz + 1;
            if (nz>=0 && nz<NZspl){
                
                betaZ = betaZ_pc[lz] ;
                for (lx = 0 ; lx<4 ; lx++  ){
                    
                    nx = nghbCorner[1] + lx + 1;
                    if (nx>=0 && nx<NXspl){
                        
                        betaX = betaX_pc[lx] ;
                        for (ly = 0 ; ly<4 ; ly++  ){
                            ny = nghbCorner[0] + ly + 1;
                            
                            if (ny>=0 && ny<NYspl){
                                betaY = betaY_pc[ly] ;
                                
                                lll = ny + nx*NYspl +  nz*NXspl*NYspl ;
                                bbb = betaX * betaY * betaZ ;
                                
                                
                                if (bbb>0){ // this is to guarantee the matrix do not contain superfluous zeros 
                                    //printf("bbb= %f\n",bbb) ;
                                    mxSplineGrad[entry_counter] = bbb  ;
                                    ir[entry_counter] = (mwIndex)(lll) ;
                                    entry_counter++ ;
                                }
                            
                            }
                            
                        }
                        
                    }
                    
                }
                
            }
            
        }

    }
    
    jc[   (mwIndex)(Ngrid_all)   ] =  entry_counter  ; // see Matlab documentation... 
}













