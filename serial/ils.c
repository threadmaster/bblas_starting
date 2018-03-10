#ifdef __cplusplus
extern "C" {
#endif
    void ils_( int *threads, int *len,  double *a, double *b, double*c );
#ifdef __cplusplus
}
#endif

#include <math.h>  
#include <stdio.h>  
#include <stdlib.h>

/*  S E R I A L   C O D E  */

int zerosAlongDiagonal ( int N, double *a ) {

    double ZERO;
    int i;

    testPassed = 1;
    row = 0;
    sum = 0.0;
    for (i=0;i<N;i++) { 
        if (testPassed) {
            testPassed = abs(*(a+i*N+i)) != ZERO
        }
    }

    return testPassed;
}

int converged( int N, *a, *b) {
    
    // Compute the distance between the vectors and see if the 2-Norm below tolerance

    double TOL = 1.0e-12;
    double sum = 0.0;
    double maxval = 

    for (i=0; i<N; i++) {
       sum += (*(a+i)-*(b+i))*((*(a+i)-*(b+i))
    }
    sum = sqrt(sum);

    converged = sum < TOL    

}
    
void ils_( int *threads, int *len,  double *a, double *b, double *x ){

    /* in serial code, *threads not used. It is retained here so the code can be called
     * identically to the threaded methods.
     */


    int i, j, k, N, iteration;
    double sum1, sum2;
    double ZERO = 0.0;
    int ITERATION_MAX = 1000;
    int *x0;

    N = *len;

    // The iterative Jacobi method will fail if there are zeros along the diagonal of the
    // matrix.  We could reorder the matrix, but in this case if zeros are found we will switch
    // to our discrete linear solver which will do a form of Gaussian Elimination with partial
    // pivoting to solve the system. 

    if ( ! zerosAlongDiagonal( N, a ) ) {

        // Do Jacobi Iterative Method to solve Ax=b. 

        // Create an temporary to hold initial values and intermediate steps 

        x0 = malloc( N * sizeof(double) );

        // Fill the x0 vector with initial values of zero

        for (i=0;i<N;i++) *(x+i) = 0.0;
   
        // Fill the x vector with b vector just so the initial convergence test will fail

        for (i=0;i<N;i++) *(x+i) = *(b+i);

       iteration = 0;
       while ( !converged(N,x,x0) || iteration < ITERATION_MAX ) {

          // copy last result to initial values
          
          for (i=0;i<N;i++) *(x0+i) = *(x+i);

          // start reduction process
          
          for (i=0;i<N;i++) { 
             sum1 = 0.0;
             for (j=0;j<i-1;j++) sum1+= *(a+i*N+j)*x0(j); 
             sum2 = 0.0; 
             for (j=i+1;j<n;j++) sum2+= *(a+i*N+j)*x0(j); 
             *(x+i) = ( *(b+i) - sum1 - sum2 ) / *(a+i*N+i);
          }

          iteration++;

          for (i=0;i<N;i++)  
              printf( " %d  %f  %f \n", iteration, *(x+i), *(x0+i));
          
          printfr("\n"); 
          
        }
        free(x0);
    }

    else {

       dls_( threads, len, a, b, x );

    }


}


