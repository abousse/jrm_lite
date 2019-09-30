// JRM_Lite
// Copyright University College London 2015, 2016, 2017, 2018, 2019
// Author: Alexandre Bousse, Institute of Nuclear Medicine, UCL
// For research purpose only.

#include "jrm_source_mex.h"

// mex  -largeArrayDims pm_warp3Dderivatives_mex.cpp

// Jacobian of Wu

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
    // -----------------------------------------------------------------------------
    // INPUTS
    // -----------------------------------------------------------------------------
    
    
    double *mxVol = mxGetPr(prhs[0]) ; // volume to warp
    
    double *mxXnew = mxGetPr(prhs[1]) ;  // x-coordinates in the original grid
    double *mxYnew = mxGetPr(prhs[2]) ;  // y-coordinates in the original grid
    double *mxZnew = mxGetPr(prhs[3]) ;  // z-coordinates in the original grid
    
    // -----------------------------------------------------------------------------
    
    mwSize const *dimsGrid = mxGetDimensions(prhs[1]) ; // image dimensions 
    mwSize const NX = dimsGrid[1] ; 
    mwSize const NY = dimsGrid[0] ; 
    mwSize const NZ = dimsGrid[2] ; 
    mwSize const Ngrid_all = NX*NY*NZ ;// length of mxXgrid, i.e. number of voxels
    
    
    int j, l,ll, k, lll ;
    double Xnew, Ynew, Znew, Xdiff, Ydiff, Zdiff ;
    
    int lz, lx, ly ;
    int nz, nx, ny ;
    double betaX, betaY, betaZ  ;
    double betaDerX, betaDerY, betaDerZ  ;
    
    double betaX_pc[4], betaY_pc[4], betaZ_pc[4] ; // pre-computed values
    double betaDerX_pc[4], betaDerY_pc[4], betaDerZ_pc[4] ;
    int nghbCorner[3] ;
    
    
    double Dbb, bDb, bbD, Dww, wDw, wwD, val ; // used for beta'*beta*beta,  etc. (see paper)
    
    
    //-----------------------------------------------------------------------------------------------------------------------
    // output: mxDww, mxwDw, mxwwD, ie diag[W_X (alpha)u],   diag[W_Y (alpha)u] and  diag[W_Z (alpha)u]
    //-----------------------------------------------------------------------------------------------------------------------
    
    double *mxDww = NULL  ;
    double *mxwDw = NULL  ;
    double *mxwwD = NULL  ;
    
    mwIndex *ir0 = NULL ;
    mwIndex *jc0 = NULL ;
    mwIndex *ir1 = NULL ;
    mwIndex *jc1 = NULL ;
    mwIndex *ir2 = NULL ;
    mwIndex *jc2 = NULL ;
    
    mwIndex entry_counter = 0 ;
    
    plhs[0] = mxCreateSparse(Ngrid_all,Ngrid_all,Ngrid_all, mxREAL) ; mxDww = mxGetPr(plhs[0]) ; // D = derivative
    plhs[1] = mxCreateSparse(Ngrid_all,Ngrid_all,Ngrid_all, mxREAL) ; mxwDw = mxGetPr(plhs[1]) ;
    plhs[2] = mxCreateSparse(Ngrid_all,Ngrid_all,Ngrid_all, mxREAL) ; mxwwD = mxGetPr(plhs[2]) ;
    
    ir0 = mxGetIr(plhs[0]) ;
    jc0 = mxGetJc(plhs[0]) ;
    
    ir1 = mxGetIr(plhs[1]) ;
    jc1 = mxGetJc(plhs[1]) ;
    
    ir2 = mxGetIr(plhs[2]) ;
    jc2 = mxGetJc(plhs[2]) ;
    
    
    /////////////////////////////////////////////////////////////
    
    for (j=0; j<Ngrid_all; j++){  // can be parallelised 
        
        Dww = 0 ; wDw = 0 ; wwD = 0 ;
        
        jc0[(mwIndex)(j)] = entry_counter ;
        jc1[(mwIndex)(j)] = entry_counter ;
        jc2[(mwIndex)(j)] = entry_counter ;
        
        Xnew = mxXnew[j] ; Ynew = mxYnew[j] ; Znew = mxZnew[j] ;
        
        
        initInt(nghbCorner,3) ;
        initDouble(betaX_pc,4) ; initDouble(betaY_pc,4) ; initDouble(betaZ_pc,4) ;
        initDouble(betaDerX_pc,4) ; initDouble(betaDerY_pc,4) ; initDouble(betaDerZ_pc,4) ;
        findCorner(Xnew, Ynew, Znew, 1,1,1, nghbCorner) ;
        
        for (ll = 0; ll<4 ; ll++ ){
            nz = nghbCorner[2] + ll ;
            nx = nghbCorner[1] + ll ;
            ny = nghbCorner[0] + ll ;
            
            Zdiff = Znew - nz ;
            betaZ = beta(Zdiff) ; betaDerZ = betaDer(Zdiff) ;
            betaZ_pc[ll] = betaZ ; betaDerZ_pc[ll] = betaDerZ ;
            Xdiff = Xnew - nx ;
            betaX = beta(Xdiff) ; betaDerX = betaDer(Xdiff) ;
            betaX_pc[ll] = betaX ; betaDerX_pc[ll] = betaDerX ;
            Ydiff = Ynew - ny ;
            betaY = beta(Ydiff) ; betaDerY = betaDer(Ydiff) ;
            betaY_pc[ll] = betaY ; betaDerY_pc[ll] = betaDerY ;
        }
        
        
        
        for (lz = 0; lz<4 ; lz++){
            
            nz = nghbCorner[2] + lz ;
            if (nz>=0 && nz<NZ){
                
                betaZ = betaZ_pc[lz] ;
                betaDerZ = betaDerZ_pc[lz] ;
                
                for (lx = 0 ; lx<4 ; lx++  ){
                    
                    nx = nghbCorner[1] + lx ;
                    if (nx>=0 && nx<NX){
                        
                        betaX = betaX_pc[lx] ;
                        betaDerX = betaDerX_pc[lx] ;
                        
                        for (ly = 0 ; ly<4 ; ly++  ){
                            ny = nghbCorner[0] + ly ;
                            
                            if (ny>=0 && ny<NY){
                                
                                betaY = betaY_pc[ly] ;
                                betaDerY = betaDerY_pc[ly] ;
                                
                                lll = ny + nx*NY +  nz*NX*NY ;
                                val = mxVol[lll] ;
                                
                                Dbb = betaDerX * betaY * betaZ ;
                                bDb = betaX * betaDerY * betaZ ;
                                bbD = betaX * betaY * betaDerZ ;
                                
                                Dww = Dww + Dbb*val ;
                                wDw = wDw + bDb*val ;
                                wwD = wwD + bbD*val ;
                                
                            }
                            
                        }
                        
                    }
                    
                }
                
            }
            
        }
        // /!\ WARNING /!\
        // it seems that some of the entries are still zeros, despite the fact that the sum is performed over neighbours 
        // that should be in the support of the Bspline. Maybe it's because of vol[lll]?
        // edit - it is likely to be because of vol[lll]. This can be solved by multiplying the sparse output
        // by a scalar. However, a neat implementation requieres  3 counters for Dww, wDw and wwD and test if each value is non-zero
        
        mxDww[entry_counter] = Dww ;
        mxwDw[entry_counter] = wDw ;
        mxwwD[entry_counter] = wwD ;
        
        
        ir0[entry_counter] = (mwIndex)(j) ;
        ir1[entry_counter] = (mwIndex)(j) ;
        ir2[entry_counter] = (mwIndex)(j) ;
        
        entry_counter++ ;
    }
    
    jc0[   (mwIndex)(Ngrid_all)   ] =  entry_counter  ;
    jc1[   (mwIndex)(Ngrid_all)   ] =  entry_counter  ;
    jc2[   (mwIndex)(Ngrid_all)   ] =  entry_counter  ;
    
}



