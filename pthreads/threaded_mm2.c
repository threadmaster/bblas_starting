#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#define ROWS 1000 
#define COLUMNS 1000 
#define NTHREADS ROWS*COLUMNS 


void *thread_function();
void mmm( int, int, *double, *double, *double );

double *A, *B, *C;


struct args {
    int N; 
    int row;
    int column;
    double *Aptr;
    double *Bptr;
    double *Cptr;
};


 
main()
{
   pthread_t *thread_id;
   struct args *thread_args;
   int i, j;
   double trace;

   A =  ( double * ) malloc( NTHREADS * sizeof(double) );
   B =  ( double * ) malloc( NTHREADS * sizeof(double) );
   C =  ( double * ) malloc( NTHREADS * sizeof(double) );

   thread_id = (pthread_t *) malloc (NTHREADS * sizeof(pthread_t));

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
  
  
   
   for(i=0; i < ROWS; i++) {
      for(j=0; j < COLUMNS; j++)
      {
         thread_args = ( struct args * )  malloc(sizeof( struct args));
//         printf("---->  %i  %i  %i  %i\n", NTHREADS, ROWS, i, j); 
         thread_args->N   = ROWS;
         thread_args->row = i;
         thread_args->column = j; 
         thread_args->Aptr = A;
         thread_args->Bptr = B;
         thread_args->Cptr = C;

         pthread_create( thread_id+(ROWS*i+j), NULL, &thread_function, thread_args );
      }
   }

   
   for (i=0; i < ROWS; i++) {
      for(j=0; j < COLUMNS; j++)
      {
         pthread_join( thread_id[ROWS*i+j], NULL); 
      }
   }
  
   /* Now that all threads are complete I can print the final result.     */
   /* Without the join I could be printing a value before all the threads */
   /* have been completed.                                                */

  
   trace = 0.0;
   for(i=0; i < ROWS; i++) {
      for(j=0; j < COLUMNS; j++)
      {
         if ( i == j ) { trace = trace + *(C+ROWS*i+j);}
         //printf(" %f ", *(C+(ROWS*i+j)));
      }
      //printf("\n", NULL );
   }
      printf("\n\n Trace = %f\n\n ", trace );
}

void *thread_function( struct args *thread_args  )
{
   int k;
   double val;
   int row, column, veclen;
   double *AMatrix, *BMatrix, *CMatrix;

   veclen = thread_args->N;
   row    = thread_args->row;
   column = thread_args->column;
   AMatrix  = thread_args->Aptr;
   BMatrix  = thread_args->Bptr;
   CMatrix  = thread_args->Cptr;

   val = 0.0;
   for (k=0;k<veclen;k++) {
      val = val + *(AMatrix+(row*veclen+k)) * *(BMatrix+(column*veclen+k));
   }
//   printf(" In Function --  %i  %i   %i   %f\n", veclen, row, column, val);

   *(CMatrix+(row*veclen+column)) = val;

   free(thread_args);
   pthread_exit(NULL);
}

   void mmm( int numThreads, int matrixDimension, *double A, *double B, *double C ){

   // This progran has to break up the data, spawn the processes, gather the results, and 
   // clean up after itself.

   // if thera are fewer dimensions than threads, do a simple matrix multiplication.

   if ( matrixDimension < numThreads ) {
     printf("Do normal matrix multipy here.");
   else

      
