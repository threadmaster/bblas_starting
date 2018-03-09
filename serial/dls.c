#ifdef __cplusplus
extern "C" {
#endif
    void mmm_( int *threads, int *len,  double *a, double *b, double*c );
#ifdef __cplusplus
    }
#endif

#include <math.h>  
#include <stdio.h>  
#include <stdlib.h>

/*  S E R I A L   C O D E  */

int strictlyDiagonallyDominant( int N, double *a ) {

    double sum;
    int i, testPassed, row;

    testPassed = 1;
    row = 0;
    sum = 0.0;
    for (row=0;row<N;row++) { 
        if (testPassed) {
            sum = 0.0;
            for (i=0;i<N;i++) sum+=*(a+row*N+i);
            sum-=fabs(*(a+row*N+row)); 
            testPassed = fabs(*(a+row*N+row)) > sum;
        }
    }

    printf("-- diagonal dominance test completed, status = %d\n", testPassed);

    return testPassed;
}
        

void dls_( int *threads, int *len,  double *a, double *b, double *x ){

/* in serial code, *threads not used. It is retained here so the code can be called
 * identically to the threaded methods.
 */


    int i, j, k, N, u;
    int singular, iPivot, rows, rows2;
    double pivotMax, tmp, *y;
    double sum;
    double ZERO = 0.0;
    int *p;
   
    N = *len;
 
// Check A for strict diagonal dominance to see if we can reduce the matrix without 
// doing any row interchanges.   We could also check for positive definiteness to
// achieve the same thing.

    if ( ! strictlyDiagonallyDominant( N, a ) ) {

        // Do Gaussian Elimination with Partial Pivoting 

        // Create an array to hold pivot swaps 

        p = malloc( (N-1) * sizeof(int) );

        for (k=0;k<N-1;k++) *(p+k)=k;

        // Search for largest value in the column and swap the 
        // entire row containing that value with the current
        // pivot row.

        printf("Starting A matrix...\n"); 
        // print A
        for (i=0;i<N;i++) {
            for (j=0;j<N;j++) {
                printf(" %10.5f ", *(a+N*i+j));
            }
            printf("\n");
        }

        for (k=0;k<N-1;k++) {
            pivotMax = *(a+k*N+k);
            iPivot = k; 
            for (u=k;u<N;u++) {
                if ( fabs(*(a+u*N+k)) > fabs(pivotMax) ) {
                    pivotMax = *(a+u*N+k);
                    iPivot = u;
                }
            }
            // If a greater pivot value was found, swap the rows.
            if ( iPivot != k ) {
                printf("Swapping rows %d and %d\n", k, iPivot);
                u = iPivot; 
                for (j=k;j<N;j++) {
                    tmp = *(a+k*N+j);
                    *(a+k*N+j) = *(a+u*N+j);
                    *(a+u*N+j)=tmp;
                }
            }
            *(p+k) = iPivot;
            if ( *(a+k*N+k) != ZERO ) {
                for (rows=k+1;rows<N;rows++) { 
                    *(a+rows*N+k) = *(a+rows*N+k) / *(a+k*N+k);

                for (rows2=k+1;rows2<N;rows2++) { 
                    *(a+rows*N+rows2) = *(a+rows*N+rows2) - 
                        *(a+rows*N+k) * *(a+k*N+rows2) ;
                }
               }

                    for (i=0;i<N;i++) {
                        for (j=0;j<N;j++) {
                            printf("++++ %10.5f ", *(a+N*i+j));
                        }
                        printf("\n");
                    }
                    printf("\n");
                }
            else {

                /* Handle the case of a zero pivot element, singular matrix */

                printf( " *** MATRIX A IS SINGULAR *** \n");
                printf( "    -- EXECUTION HALTED --\n");
                exit(1);
            }

            }
            // Now that we know we have reduced the matrices, start back
            // substitution process to solve for vector x.

            // print A
            for (i=0;i<N;i++) {
                for (j=0;j<N;j++) {
                    printf(" %10.5f ", *(a+N*i+j));
                }
                printf("\n");
            }
            printf("\n");

        /* We now need to arrange b so that it has undergone the 
         * operations as the matrix a.  This will then form
         * the vector y for the solution of Ux=y where U is the 
         * upper-triangular matrix formed in the elimination process
         * above. 
         */

        for (k=0; k<N-1; k++ ) {
            // Swap rows x with p(k) 
            tmp = *(b+k);
            *(b+k) = *(b+ *(p+k));
            *(b+ *(p+k)) = tmp;
            
            for (j=k+1;j<N;j++) 
                *(b+j)= *(b+j) - *(b+k) * *(a+N*j+k);  
        } 

        // Now we can do the backward substitution to get the solution
        // vector x
     
        *(b+N-1) = *(b+N-1) / *(a+N*(N-1)+(N-1));
        for (i=N-2;i>=0;i--){
           tmp = 0.0;
           for (j=i+1;j<N;j++) {
             tmp = tmp + *(a+i*N+j) * *(b+j);
           }
           *(b+i) = ( *(b+i) - tmp ) / *(a+i*N+i); 
        }

        for (i=0;i<N;i++) *(x+i) = *(b+i);
        printf("Solutions vector\n");
        for (k=0; k<N; k++ ) {
            printf( "x[%d] = %f\n", k, *(x+k));
        } 
    }

           
    else {

        // Since we know the matrix is diagonally dominant, verify
        // that none of the pivot elements are equal to zero

        singular = 1; 
        while ( i<N  && singular ) {
            singular = *(a+i*N+i) == ZERO;   
            i++;
        }

        if ( singular ) {
            printf( " *** MATRIX A IS SINGULAR *** \n");
            printf( "    -- EXECUTION HALTED -- \n");
            exit(1);
        }

        // Now we know that we have a diagonally dominant matrix that is
        // not singular -- so it sould be possible to do the LU factorization
        

        for (k=0; k<N-1; k++) {
            for (rows=k+1;rows<N;rows++) {
                *(a+rows*N+k) = *(a+rows*N+k) / *(a+k*N+k);

                for (rows2=k+1;rows2<N;rows2++) { 
                    *(a+rows*N+rows2) = *(a+rows*N+rows2) - 
                        *(a+rows*N+k) * *(a+k*N+rows2) ;
                }
            }
        }
        
        // At this point the LU factorizaton should be done and we have to do two
        // triangular back substitutions.  The solution to Ax=b is solved by first 
        // solving Ly=b for y and then Ux=y for the solution vector x.

        // Solving lower-triangular (Ly=b) first

        for (k=0; k<N-1; k++ ) {
            for (j=k+1;j<N;j++) 
                *(b+j)= *(b+j) - *(b+k) * *(a+N*j+k);  
        } 

        // Now we can do the backward substitution to get the solution
        // vector x for the upper-triangular system (Ux=y)
     
        *(b+N-1) = *(b+N-1) / *(a+N*(N-1)+(N-1));
        for (i=N-2;i>=0;i--){
           tmp = 0.0;
           for (j=i+1;j<N;j++) {
             tmp = tmp + *(a+i*N+j) * *(b+j);
           }
           *(b+i) = ( *(b+i) - tmp ) / *(a+i*N+i); 
        }

        for (i=0;i<N;i++) *(x+i) = *(b+i);

        // At this point the solution to the system should be in vector x 


    }

}


