#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#define ROWS 4000 
#define COLUMNS 4000 
#define NTHREADS ROWS*COLUMNS 


void *mmm_thread_worker();

double *A, *B, *C;


struct args {
    int N; 
    int startRow;
    int stopRow;
    double *Aptr;
    double *Bptr;
    double *Cptr;
};

void mmm( int numThreads, int matrixDimension, double *A, double *B, double *C ){

    // This progran has to break up the data, spawn the processes, gather the results, and 
    // clean up after itself.


    int *numberOfRows;
    int startRow, stopRow;
    pthread_t *thread_id;
    struct args *thread_args;


    // if there are fewer dimensions than threads, do a simple matrix multiplication.
    if ( matrixDimension < numThreads ) {
       for (int i=0;i<matrixDimension;i++) {
           for (int j=0;j<matrixDimension;j++) {
               *(C+(i*matrixDimension+j))=0.0;
               for (int  k=0;k<matrixDimension;k++) {
                   *(C+(i*matrixDimension+j)) += *(A+(i*matrixDimension+k)) * *(B+(k*matrixDimension+j));
               }
           }  
       } 
    }

    else { 

    /* 
     * Now the parallel work begins.  The process is to first determine how to break
     *  up the matrix work equitably across the threads.  Once this is done the struct is filled with
     *  the information and a thread is started using the information.  Other than the size of the
     *  matrices and the rows to be processed, only the pointers to the locations in memory of the matrices
     *  are passed.
     */ 

        // Malloc an array to keep up with thread id's for each thread
        thread_id = (pthread_t *) malloc (NTHREADS * sizeof(pthread_t));

        // Malloc an array to keep up with how many rows to work on in each thread
        numberOfRows = ( int * ) malloc( NTHREADS * sizeof(int) );

        for (int i=0; i<numThreads; i++ ){
            *(numberOfRows+i) = matrixDimension / numThreads;
        }
        for (int i=0; i< matrixDimension % numThreads; i++ ){
            *(numberOfRows+i) = *(numberOfRows+i) + 1;
        }

        //        for (int i=0; i<numThreads; i++ ){
        //           printf("Thread %d will handle %d rows\n", i, *(numberOfRows+i) );

        stopRow=0;
        for(int i=0; i < numThreads ; i++) {
            {   
                startRow=stopRow;
                stopRow=startRow+*(numberOfRows+i);
                thread_args = ( struct args * )  malloc(sizeof( struct args));
                thread_args->N   = matrixDimension;
                thread_args->startRow = startRow;
                thread_args->stopRow = stopRow; 
                thread_args->Aptr = A;
                thread_args->Bptr = B;
                thread_args->Cptr = C;

                pthread_create( thread_id+i, NULL, &mmm_thread_worker, thread_args );
            }
        }
        for(int i=0; i < numThreads ; i++) {
          pthread_join( *(thread_id+i), NULL); 
        }

        free(numberOfRows);
        free(thread_id);
    }

    }


    int main()
    {
    //    pthread_t *thread_id;
    //    struct args *thread_args;
        int i, j;
        double trace;

        A =  ( double * ) malloc( NTHREADS * sizeof(double) );
        B =  ( double * ) malloc( NTHREADS * sizeof(double) );
        C =  ( double * ) malloc( NTHREADS * sizeof(double) );

     //   thread_id = (pthread_t *) malloc (NTHREADS * sizeof(pthread_t));

        /* Fill A and B matrices */

        for(i=0; i < ROWS; i++) {
            for(j=0; j < COLUMNS; j++) {
                if (i != j)  {
                    *(A+(ROWS*i+j)) = 0.0; 
                    *(B+(ROWS*i+j)) = 0.0; 
                }
                else {
                    *(A+(ROWS*i+j)) = 1.0; 
                    *(B+(ROWS*i+j)) = 1.0; 
                }
            }
        }

        mmm( 8, ROWS, A, B, C );  

 
        printf("Matrix multiply should be done.\n");

        if (ROWS < 8 ) {
        for (int i=0; i<ROWS; i++) {
          for (int j=0; j<ROWS; j++ ) 
            printf(" %f ", *(C+i*ROWS+j) );
          printf("\n");
        }
        }

        // compute sum of diagonal elements
        trace = 0.0;
        for (i=0; i<ROWS; i++) trace += *(C+(i*ROWS+i));

        printf("The dimension is %d and the trace is %f\n", ROWS, trace);

}        

       
       void *mmm_thread_worker( struct args *thread_args  ) {

       int i, j, k;
       double val;
       int rowStart, rowStop, N; 
       double *A, *B, *C;
        
       N        =  thread_args->N;
       rowStart =  thread_args->startRow;
       rowStop  =  thread_args->stopRow; 
       A        =  thread_args->Aptr;
       B        =  thread_args->Bptr;
       C        =  thread_args->Cptr;

       for (i=rowStart;i<rowStop;i++) {
           for (j=0;j<N;j++) {
               *(C+(i*N+j))=0.0;
               for (k=0;k<N;k++) {
                   *(C+(i*N+j)) += *(A+(i*N+k)) * *(B+(k*N+j));
               }
           }  
       } 

     free(thread_args);
     pthread_exit(NULL);
     }
     



