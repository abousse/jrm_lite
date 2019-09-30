// JRM_Lite
// Copyright University College London 2015, 2016, 2017, 2018, 2019
// Author: Alexandre Bousse, Institute of Nuclear Medicine, UCL
// For research purpose only.


#include "mex.h"
#include <algorithm>
#include "math.h"

struct sparseMat{
    
    mwIndex nEntries ;
    double *val ;
    mwIndex *index;
    mwSize nRow ;
    mwSize nCol ;
    
} ;

inline void initInt(int *array, int length) ;
inline void initDouble(double *array, int length) ;
inline double beta (double t) ;
inline double betaDer (double t) ;
inline void findCorner(double X, double Y, double Z, double gridWidthX, double gridWidthY, double gridWidthZ, int *nghbCorner) ;
inline int findNeighbours(double X, double Y, double Z, double gridWidth, int Ngrid, int *NghbX, int *NghbY, int *NghbZ ) ;
inline sparseMat makeJacBspline(double *Xgrid, double *Ygrid, double *Zgrid, 
        double *Xspl, double *Yspl, double *Zspl, mwSize Ngrid, mwSize Nspl, 
        double *val, mwSize *index, mwSize nEntriesMax)   ;

#include "jrm_source_mex.cpp" 



