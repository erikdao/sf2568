C Reaction-diffusion equation in 1D
C Solution by Jacobi iteration
C simple MPI implementation
C
C C Michael Hanke 2006-12-14
C
      PROGRAM poisson1D
      IMPLICIT NONE
      INTEGER nn, smx
C
C define problem to be solved
C   nn:  number of inner grid points
C   smx: number of iterations
      PARAMETER(nn = 1000, smx = 1000000)
C
C Use MPI
      INCLUDE 'mpif.h'
      INTEGER status(MPI_STATUS_SIZE), mpierr, signal
C
C external functions are used for r and f
      REAL*8 r, f
      EXTERNAL r, f
C
C We assume linear data distribution. The formulae according to the lecture
C  are:
C     L = N/P; (N is nn, P is pp)
C     R = N%P;
C     I = (N+P-p-1)/P;    (number of local elements) (I is ii)
C     n = p*L+MIN(p,R)+i; (global index for given (p,i)
C  Attention: We use a small trick for introducing the boundary conditions:
C     - The first ghost point on p = 0 holds u(0);
C     - the last ghost point on p = P-1 holds u(1).
C  Hence, all local vectors hold I elements while u has I+2 elements.
C
C the following arrays are much overdimensioned
      REAL*8 u(nn), unew(nn), rr(nn), ff(nn)
C
C ********** executable statements
C
C Initialize MPI
      CALL MPI_Init(mpierr)
      CALL MPI_Comm_size(MPI_COMM_WORLD, pp, mpierr)
      CALL MPI_Comm_rank(MPI_COMM_WORLD, p, mpierr)
      IF (nn .LT. pp) THEN
         WRITE(*,*) 'Too few discretization points...'
         STOP 1
      ENDIF
C
C Compute local indices for data distribution
C
C set initial guess
      DO i = 1,ii+2
         u(i) = 0d0
      ENDDO
C
C Jacobi iteration
      DO step = 1,smx
C
C RB communication
C
C local iteration step
C             unew(i) = (u(i)+u(i+2)-h*h*ff(i))/(2.0-h*h*rr(i))
       ENDDO
C
C output for graphical representation */
C Instead of using gather (which may lead to excessive memory requirements
C  on the master process) each process will write its own data portion. This
C  introduces a sequentialization: the hard disk can only write (efficiently)
C  sequentially. Therefore, we use the following strategy:
C  1. The master process writes its portion.
C  2. The master sends a signal to process 1 to start writing.
C  3. Process p waites for the signal from process p-1 to arrive.
C  4. Process p writes its portion to disk.
C  5. process p sends the signal to process p+1 (if it exists).
C
C *******************
C
C That's it
       CALL MPI_Finalize(mpierr);
       END
