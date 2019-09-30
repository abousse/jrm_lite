// JRM_Lite
// Copyright University College London 2015, 2016, 2017, 2018, 2019
// Author: Alexandre Bousse, Institute of Nuclear Medicine, UCL
// For research purpose only.


#include "jrm_source_mex.h"



void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
    double *mxVol = mxGetPr(prhs[0]) ;
    double mxWspl = mxGetScalar(prhs[1]) ;
    
    
    mwSize const *dimsVol = mxGetDimensions(prhs[0]) ; // image dimensions
    
    mwSize const N1 = dimsVol[0] ;
    mwSize const N2 = dimsVol[1] ;
    mwSize const N3 = dimsVol[2] ;
    
    
    
    
    
    int i,j,k,l ;
    int leftNeighbour ;
    int neighbours[4] ;
    double beta_pc[4] ;
    double(x) ;
    int indexVol, indexVolBinned ;
    mwSize const N1new = (int)floor(double(N1)/mxWspl) + 1 ;
    
    // double const factor = (double)N1 / (double)N1new ;
    
    //---------------------------------------------------
    // output:
    //---------------------------------------------------
    
    mwSize const dimsBinning[3] = {N1new,N2,N3} ;
    double *mxVolBinned= NULL  ;
    
    plhs[0] = mxCreateNumericArray(3, dimsBinning, mxDOUBLE_CLASS, mxREAL) ;
    mxVolBinned = mxGetPr(plhs[0]) ;
    initDouble(mxVolBinned,N1new*N2*N3)  ;
    //-----------------------------------------------------
    
    for (i=0 ; i<N1 ;i++){
        //printf("i= %i\n",i) ;
        x = i/mxWspl ;
        leftNeighbour = floor(x) - 1 ;
        
        neighbours[0] = leftNeighbour ;
        neighbours[1] = leftNeighbour + 1 ;
        neighbours[2] = leftNeighbour + 2 ;
        neighbours[3] = leftNeighbour + 3 ;
        
        
        for (l = 0; l<4 ; l++ ){
            beta_pc[l] = beta(leftNeighbour + l - x) ;
        }
        
        for (k=0 ; k<N3 ; k++){
            
            for (j=0 ; j<N2 ; j++){
                
                indexVol = i + j*N1 +  k*N1*N2 ;
                //printf("indexVol= %i\n",indexVol) ;
                
                for (l = 0; l<4 ; l++ ){
                    
                    if (neighbours[l]>=0 && neighbours[l]<N1new){
                        
                        indexVolBinned = neighbours[l] + j*N1new +  k*N1new*N2 ;
                        //printf("indexVolBinned= %i\n",indexVolBinned) ;
                        mxVolBinned[indexVolBinned] = mxVolBinned[indexVolBinned] + mxVol[indexVol]*beta_pc[l] ;
                        
                    }
                }
                
            }
            
        }
        
    }
    
}
