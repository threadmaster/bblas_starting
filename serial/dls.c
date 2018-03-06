#ifdef __cplusplus
extern "C" {
#endif
    void mmm_( int *threads, int *len,  double *a, double *b, double*c );
#ifdef __cplusplus
    }
#endif

/*  S E R I A L   C O D E  */


int strictlyDiagonallyDominant( int N, double *a ) {

    double sum;
    int i, testPassed, row, N;

    N = *len;

    testPassed = 1;
    row = 0;
    sum = 0.0;
    for (row=0;row<N;i++) { 
        if (testPassed) {
            sum = 0.0;
            for (i=0;i<N;i++) sum+=*(a+row*N+i);
            sum-=fabs(*(a+row*N+row)); 
            testPassed = fabs(*(a+row*N+row)) > sum;
        }
    }

    return testPassed;
}
        

void dls_( int *threads, int *len,  double *a, double *b, double *x ){

/* in serial code, *threads not used. It is retained here so the code can be called
 * identically to the threaded methods.
 */


    int i, j, k, N, u;
    int singular;
    double pivotMax, tmp, *y;
    double ZERO = 0.0;
   
    N = *len;
 
// Check A for strict diagonal dominance to see if we can reduce the matrix without 
// doing any row interchanges.   We could also check for positive definiteness to
// achieve the same thing.

    if ( ! strictlyDiagonallyDominant( N, a ) ) {
 
        // Do Gaussian Elimination with Partial Pivoting 
        
       // Search for largest value in the column and swap the 
       // entire row containing that value with the current
       // pivot row.
       for (k=0;k<N;k++) {
          pivotMax = *(a+i*N+i);
          iPivot = i 
          for (k=i;k<N;k++) {
              if ( fabs(*(a+k*N+i)) > fabs(pivot) ) {
                 pivot = *(a+k*N+i);
                 iPivot = k;
              }
          }
       // If a greater pivot value was found, swap the rows.
          if ( iPivot != i ) {
              *(swaps+i) -> from = iPivot;
              *(swaps+i) -> to   = i;
              for (k=0;k<N;k++) {
                 tmp = *(a+i*N+k);
                 *(a+i*N+k) = *(a+iPivot*N+k);
                 *(a+iPivot*N+k)=tmp;
              }
          }
       }

       
    }       
      
    
    else {

       // Since we know the matrix is diagonally dominant, verify
       // that none of the pivot elements is equal to zero

       singular = 1 
       while ( i<N  && singular ) {
          singular = *(a+i*N+i) == ZERO;   
          i++;
       }
    
       if ( singular ) {
         puts( " *** MATRIX A IS SINGULAR *** ");
         puts( "    -- EXECUTION HALTED --");
         exit(1);
       }

       // Now we know we that have a diagonally dominant matrix that is
       // not singular -- so it sould be possible to do the LU factorization

       for (k=0; k<N-1; k++) {
          for (rows=k+1;rows<N;rows++) {
             *(a+rows*N+k) = *(a+rows*N+k) / *(a+k*rows+k);
             *(a+rows*N+rows) = *(a*rows*N+rows) - *(a+rows*N+k) * *(a+k*N+rows);
          }
       }

       // At this point the LU factorizaton should be done and we have to do two
       // triangular back substitutions.  The solution to Ax=b is solved by first 
       // solving Ly=b for y and then Ux=y for the solution vector x.

       // Solving lower-triangular (Ly=b) first

       y = malloc( N * sizeof(double) );
       
       *(y+0) = *(b+0) / *(a+0*N+0)
       for (i=1;i<N;i++) {
          sum = 0.0;
          for (j=0;j<i-1;j++) {
            sum+= *(a+i*N+j) * *(y+j);
          }
          *(y+i) = ( *(b+i) - sum ) / *(a+i*N+i);
       }

       // Now solve upper-triangular (Ux=y) using results from prior step

       *(x+(N-1)) = *(x+(N-1)) / *(a+(N-1)*N+(N-1));
       for (i=N-2;i>=0;i--} {
          sum = 0.0;
          for (j=i+1;j<N;j++) {
            sum+= *(a+i*N*j) * *(x+j);
          }
          *(x+i) = ( *(y+i) - sum ) / *(a+i*N+i);
       }
 
       // At this point the solution to the system should be in vector x 
     

    }


// Normal Matrix Multiplication
/*
    for (i=0; i<veclen; i++) {
        for (j=0; j<veclen; j++) {
            *(c+(i*veclen+j)) = 0.0;
            for (k=0;k<veclen;k++){
                *(c+(i*veclen+j)) += *(a+(i*veclen+k)) * *(b+(k*veclen+j)); 
            }
        }
    }
*/

#endif
}


