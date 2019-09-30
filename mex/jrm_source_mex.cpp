// JRM_Lite
// Copyright University College London 2015, 2016, 2017, 2018, 2019
// Author: Alexandre Bousse, Institute of Nuclear Medicine, UCL
// For research purpose only.

void initInt(int *array, int length){
    int i ;
    for (i=0 ; i<length ; i++){
        array[i] = 0 ;
        /*printf("array[i] = %i\n",array[i]) ;*/
    }
}

void initDouble(double *array, int length){
    int i ;
    for (i=0 ; i<length ; i++){
        array[i] = 0 ;
        /*printf("array[i] = %f\n",array[i]) ;*/
    }
}



double beta (double t){
    double abst, abst2, abst3, abstm2,spline ;
    abst = fabs((double)t) ;
    abstm2 = abst - 2 ;
    abst2 = abst*abst ;
    abst3 = abst2*abst ;
    
    if(abst<2){
        
        if(abst<1){
            spline = 0.5*abst3 - abst2 + 2.0/3.0 ;
            return spline ;
        }
        else{
            spline = -abstm2*abstm2*abstm2 / 6 ;
            return spline ;
        }
    }
    else{
        spline = 0 ;
        return spline ;
    }
}


double betaDer (double t){
    
    double tm2,tp2,abst,spline ;
    abst = fabs((double)t) ;
    
    if (abst<2){
        
        if (t>=1){
            tm2 = t-2 ;
            spline = -0.5*tm2*tm2 ;
            return spline ;
        }
        else if (t>=0){
            spline = 1.5*t*t -2*t ;
            return spline ;
        }
        else if (t>=-1){
            spline = -1.5*t*t - 2*t ;
            return spline ;
        }
        else {
            tp2 = t+2 ;
            spline = 0.5*tp2*tp2 ;
            return spline ;
        }
    }
    else{
        spline = 0 ;
        return spline ;
    }
}

void findCorner(double X, double Y, double Z, double gridWidthX,double gridWidthY,double gridWidthZ, int *nghbCorner){
    
    int Xgrid = (int)floor(X / gridWidthX) ;
    int Ygrid = (int)floor(Y / gridWidthY) ;
    int Zgrid = (int)floor(Z / gridWidthZ) ;
    
    nghbCorner[0] = Ygrid -1;  // can be outside of the box
    nghbCorner[1] = Xgrid -1;
    nghbCorner[2] = Zgrid -1;
    
}



int findNeighbours(double X, double Y, double Z, double gridWidth, int Ngrid, int *nghbX, int *nghbY, int *nghbZ){
    // look for neighbours of [X Y Z] in a subgrid of a Matlab meshgrid of the form meshgrid(0:N-1,0:N-1) ;
    // the distance between 2 nodes of the subgrid is 'gridWidth'
    // Ngrid = Grid size where the neighbours are picked
    int nbNeighbours, lx, ly, lz, count, nx, ny, nz ;
    
    int Xgrid = (int)floor(X / gridWidth) ;
    int Ygrid = (int)floor(Y / gridWidth) ;
    int Zgrid = (int)floor(Z / gridWidth) ;
    
    // Here we seek for neighbours of X (and Y, Z)
    // X can have at most 5 neighbours (including himself). This happens only when X coincides with one node of the grid, 
    // in which case the 2 extrem neighbours (1 and 5) are exactly at the edges of the b-pline function.
    // In a regular case, X has only 4 neighbours.
    // X is floored to an integer. This integer is the first left neighbour of X. Therefore, it is only necessary to loop
    // from -1 to +2

    
    nbNeighbours = 0 ;
    for (lz=-1 ; lz<=2 ; lz++){
        
        nz = Zgrid + lz ;   
        
        if (nz>=0 && nz<Ngrid){
            
            for (lx=-1 ; lx<=2 ; lx++){
                
                nx = Xgrid  + lx ;
                
                if (nx>=0 && nx< Ngrid){
                    
                    for (ly=-1 ; ly<=2 ; ly++){
                        
                        ny = Ygrid + ly ;
                        
                        if (ny>=0 && ny<Ngrid){
                            
                            nbNeighbours++ ;
                            nghbX[nbNeighbours-1] = nx ;
                            nghbY[nbNeighbours-1] = ny ;
                            nghbZ[nbNeighbours-1] = nz ;
                            
                        }
                    }
                }
            }
        }
    }
    
    return nbNeighbours ;
}

















