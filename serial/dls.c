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
        

void dls_( int *threads, int *len,  double *a, double *b, double *c ){

/* in serial code, *threads not used. It is retained here so the code can be called
 * identically to the threaded methods.
 */


    int i, j, k, N, iPivot;
    struct pairs{
       int from;
       int to; 
    };
    pairs *swaps;
    double pivot, tmp;

   
    N = *len;
 
// Check A for strict diagonal dominance to see if we can reduce the matrix without 
// doing any row interchanges.   We could also check for positive definiteness to
// achieve the same thing.

    if ( ! strictlyDiagonallyDominant( N, a ) ) {
 
        // Do Row Interchanges
        
       swaps = malloc(N * sizeof(pairs));

       // Search for largest value in the column and swap the 
       // entire row containing that value with the current
       // pivot row.
       for (i=0;i<N;i++) {
          pivot = *(a+i*N+i);
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

       // Do Row Reductions 
      

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


