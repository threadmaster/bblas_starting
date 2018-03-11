program driver 

integer :: NDIM

real (kind=8) :: wall_start, wall_end
real (kind=8) :: cpu_start, cpu_end
real (kind=8) :: trace


integer :: startval, stopval, stepval, nthreads
real (kind=8) :: walltime
real (kind=8) :: cputime 
external walltime, cputime

character (len=8) :: carg1, carg2, carg3, carg4

real (kind=8), dimension(:), allocatable :: veca, vecb, vecx
real (kind=8), dimension(:,:), allocatable :: matrixa, matrixb, matrixc
logical :: DIAG_DOMINANT
real (kind=8) :: residual

#ifdef ACCURACY_TEST

#ifdef LS_TEST

NDIM = 100 
nthreads = 2
print *, "Performing DLS Accuracy Test"
!This portion of code is ONLY used for verifying the accuracy of the code using
!the matrix, vector b, and solution vector x stored on the class website.

!Download the files from theochem using curl (don't store these on anvil!)
!NOTE: for strictly diagonally dominant systems append _dd to last file name, e.g. -- linsolve_a_dd.dat
call system("curl -s -o linsolve_a.dat --url http://theochem.mercer.edu/csc435/data/linsolve_a.dat")
call system("curl -s -o linsolve_b.dat --url http://theochem.mercer.edu/csc435/data/linsolve_b.dat")
call system("curl -s -o linsolve_x.dat --url http://theochem.mercer.edu/csc435/data/linsolve_x.dat")

print *, "Files loaded from theochem.mercer.edu"

allocate ( matrixa(NDIM,NDIM), stat=ierr)
allocate ( veca(NDIM), stat=ierr)
allocate ( vecb(NDIM), stat=ierr)
allocate ( vecx(NDIM), stat=ierr)

open (unit=5,file="linsolve_a.dat",status="old")
do i = 1, NDIM
  do j = 1, NDIM
     read(5,*) matrixa(j,i)
  enddo
enddo
close(5)
open (unit=5,file="linsolve_b.dat",status="old")
do i = 1, NDIM
   read(5,*) vecb(i)
enddo
close(5)
open (unit=5,file="linsolve_x.dat",status="old")
do i = 1, NDIM
   read(5,*) veca(i)
enddo
close(5)

print *, "Files read into program"

! Delete the files from disk
call system("rm linsolve_a.dat linsolve_b.dat linsolve_x.dat")

print *, "Files deleted from disk."

#else

!This portion of code is ONLY used for verifying the accuracy of the code using
!the matrix and matrix inverse stored on the class website.

!Download the files from theochem using curl (don't store these on anvil!)
call system("curl -s -o matrixa.dat --url http://theochem.mercer.edu/csc435/data/matrixa.dat")
call system("curl -s -o matrixb.dat --url http://theochem.mercer.edu/csc435/data/matrixb.dat")

NDIM = 100  ! The test files are 100x100 double precision matrix and its inverse
nthreads = 2
allocate ( matrixa(NDIM,NDIM), stat=ierr)
allocate ( matrixb(NDIM,NDIM), stat=ierr)
allocate ( matrixc(NDIM,NDIM), stat=ierr)
open (unit=5,file="matrixa.dat",status="old")
do i = 1, NDIM
  do j = 1, NDIM
     read(5,*) matrixa(j,i)
  enddo
enddo
close(5)
open (unit=5,file="matrixb.dat",status="old")
do i = 1, NDIM
  do j = 1, NDIM
     read(5,*) matrixb(j,i)
  enddo
enddo
close(5)

! Delete the files from disk
call system("rm matrixa.dat matrixb.dat")
#endif


#else

! Start the normal processing here.  Read the starting, stop, and step values
! as well as the number of threads to use.
! modified to use command line arguments

call get_command_argument(1, carg1)
call get_command_argument(2, carg2)
call get_command_argument(3, carg3)
call get_command_argument(4, carg4)

! Use Fortran internal files to convert command line arguments to ints

read (carg1,'(i8)') startval
read (carg2,'(i8)') stopval
read (carg3,'(i8)') stepval
read (carg4,'(i8)') nthreads 

! Start the outermost loop to run tests of varying matrix size
 
do iter = startval, stopval, stepval

NDIM = iter

allocate ( veca(NDIM), stat=ierr)
allocate ( vecb(NDIM), stat=ierr)
allocate ( vecx(NDIM), stat=ierr)
allocate ( matrixa(NDIM,NDIM), stat=ierr)
allocate ( matrixb(NDIM,NDIM), stat=ierr)
allocate ( matrixc(NDIM,NDIM), stat=ierr)

! Build veca and vecb which, their tensor product creates the two matrices 
! to be multiplied.

do i = 1, NDIM 
     veca(i) = 1.0
     vecb(i) = 1.0 / sqrt( dble(NDIM))
enddo

! Zero the matrices using Fortran 90 syntax.
matrixa = 0.0
matrixb = 0.0

call vvm(NDIM, veca, vecb, matrixa)
call vvm(NDIM, veca, vecb, matrixb)

! If doing ILS or DLS testing, build matrix C explicitly 
! as well as solution vector X and product vector B. You
! can also specify if you want the system to be diagonally 
! dominant.

DIAG_DOMINANT = .true.

call buildLinearSystem( NDIM, matrixc, vecb, vecx,  DIAG_DOMINANT )

#endif

wall_start = walltime()
cpu_start = cputime()

!call mmm(nthreads, NDIM, matrixa, matrixb, matrixc)
!call dls(nthreads, NDIM, matrixc, vecb, vecx)
call ils(nthreads, NDIM, matrixc, vecb, vecx)

cpu_end = cputime()
wall_end = walltime()

trace = 0.0;

!do i=1, NDIM 
!     trace = trace + matrixc(i,i)
!enddo

residual = 0.0
do i=1, NDIM
   residual = residual + abs(vecx(i)-dble(i))
enddo


! Calculate megaflops based on CPU time and Walltime

! NOTE -- these need to be replaced with PAPI calls to correctly work
! for iterative linear solve
!
! Matrix multiplication is 2*N**3 flops
! Gaussian Elimination with Partial Pivoting and LU decomposition are
! approximately (2/3)*N**3 flops

! For matrix multiplication
!mflops  = 2*dble(NDIM)**3/ (cpu_end-cpu_start) / 1.0e6
!mflops2 = 2*dble(NDIM)**3/ (wall_end-wall_start)/ 1.0e6

!print *, NDIM, trace, cpu_end-cpu_start, wall_end-wall_start,  mflops, mflops2

! For discrete linear solver
mflops  = (2.0/3.0)*dble(NDIM)**3/ (cpu_end-cpu_start) / 1.0e6
mflops2 = (2.0/3.0)*dble(NDIM)**3/ (wall_end-wall_start)/ 1.0e6
 
print *, NDIM, residual, cpu_end-cpu_start, wall_end-wall_start,  mflops, mflops2


! Free the memory that was allocated based on which version of the program was
! run.

if (allocated(matrixa)) deallocate(matrixa)
if (allocated(matrixb)) deallocate(matrixb)
if (allocated(matrixc)) deallocate(matrixc)
if (allocated(veca)) deallocate(veca)
if (allocated(vecb)) deallocate(vecb)
if (allocated(vecx)) deallocate(vecx)


enddo

end program driver 


! Subroutine to build random linear systems for Solvers to Use
subroutine buildLinearSystem( N, A, B, X,  DIAG_DOMINANT )

integer :: N;
real (kind=8), dimension(N,N) :: A 
real (kind=8), dimension(N) :: B, X
logical :: DIAG_DOMINANT
real (kind=8) :: rowsum

call init_random_seed()
call random_number(A)

do i =1, N
  X(i) = dble(i)
enddo

if (DIAG_DOMINANT ) then
    do i=1,N
   rowsum = 0.0
   do j=1, N 
      rowsum=rowsum+abs(A(j,i))
   enddo
   A(i,i) = rowsum-abs(A(i,i))+10.0
 enddo
endif

B = matmul(X,A)

end subroutine 



! This is needed for the random number generator
SUBROUTINE init_random_seed() 
  INTEGER :: i, n, clock 
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed 
 
  CALL RANDOM_SEED(size = n) 
  ALLOCATE(seed(n)) 
 
  CALL SYSTEM_CLOCK(COUNT=clock) 

  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  CALL RANDOM_SEED(PUT = seed)

  DEALLOCATE(seed)
END SUBROUTINE
 
