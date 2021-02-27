/* Reaction-diffusion equation in 1D
* Solution by Jacobi iteration
* simple MPI implementation
*
* C Michael Hanke 2006-12-12
*/
#include <stdlib.h>
#include <math.h>
#include <stdio.h>


#define MIN(a,b) ((a) < (b) ? (a) : (b))

/* Use MPI */
#include <mpi.h>

/* define problem to be solved */
#define N 1000   /* number of inner grid points */
#define SMX 10000000 /* number of iterations */

/* implement coefficient functions */
extern double r(const double x);
extern double f(const double x);

double r(const double x){
  double out = -x;
  return out;
}

double f(const double x){
  double out = 2 + x * x * (x - 1);
  return out;
}


/* We assume linear data distribution. The formulae according to the lecture
are:
L = N/P;
R = N%P;
I = (N+P-p-1)/P;    (number of local elements)
n = p*L+MIN(p,R)+i; (global index for given (p,i)
Attention: We use a small trick for introducing the boundary conditions:
- The first ghost point on p = 0 holds u(0);
- the last ghost point on p = P-1 holds u(1).
Hence, all local vectors hold I elements while u has I+2 elements.
*/

int main(int argc, char *argv[])
{
  /* local variable */
  int bflag, rflag, I, step, R, P, p, i, n, L, color, tag;

  double h, x, *ff, *rr, *unew, *u;
  MPI_Status status;

  /* Initialize MPI */
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &P);
  MPI_Comm_rank(MPI_COMM_WORLD, &p);
  if (N < P) {
    fprintf(stdout, "Too few discretization points...\n");
    exit(1);
  }
  tag = 100;

  /* Compute local indices for data distribution */
  L = 1.0*N/P; // same on all processes
  R = N%P; // same on all processes
  I = (N+P-p-1)/P; // not same on all processes
  h = 1.0/(N+1); // same on all processes
  printf("I - rank %d: %d\n",p,I);
  printf("R - rank %d: %d\n",p,R);
  printf("L - rank %d: %d\n",p,L);
  printf("h - rank %d: %f\n",p,h);

  if(p==0){
    printf("number of processes: %d\n", P);
  }

  if (p % 2 == 0){
    color = 0;
    printf("Node %d: has color%d: \n",p, color );
  }else{
    color = 1;
    printf("Node %d: has color%d: \n",p, color );
  }
  /* arrays */
  unew = (double *) malloc(I*sizeof(double));
  /* Note: The following allocation includes additionally:
  - boundary conditins are set to zero
  - the initial guess is set to zero */
  u = (double *) calloc(I+2, sizeof(double));

  rr = (double *) malloc(I*sizeof(double));
  ff = (double *) malloc(I*sizeof(double));

  /* Jacobi iteration */
  for (i = 0; i < I; i++) {
    n = p*L+MIN(p,R)+i;
    x = n*h;
    ff[i] = f(x);
    rr[i] = r(x);
  }
  bflag = 0;
  rflag = 0;
  //code runs until here
  for (step = 0; step < SMX; step++) {
    /* RB communication of overlap */
    /* Black: color = 0, begins communication to the higher processor*/
    if (color == 0) {
      if (p != P-1) {
        MPI_Send(&u[I], 1, MPI_DOUBLE, p+1, tag, MPI_COMM_WORLD);
        MPI_Recv(&u[I+1], 1, MPI_DOUBLE, p+1, tag, MPI_COMM_WORLD, &status);
      }
      if (p != 0) {
        MPI_Send(&u[1], 1, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD);
        MPI_Recv(&u[0], 1, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD, &status);
      }
    }else{
      MPI_Recv(&u[0], 1, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD, &status);
      MPI_Send(&u[1], 1, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD);
      if (p != P-1) {
        MPI_Recv(&u[I+1], 1, MPI_DOUBLE, p+1, tag, MPI_COMM_WORLD, &status);
        MPI_Send(&u[I], 1, MPI_DOUBLE, p+1, tag, MPI_COMM_WORLD);
      }
    }

    for (i = 0; i < I; i++) {
      n = p*L+MIN(p,R)+i;
      x = n*h;
      /* local iteration step */
      unew[i] = (u[i]+u[i+2]-h*h*ff[i])/(2.0-h*h*rr[i]);
    }
    for (i = 0; i < I; i++) {
    u[i+1] = unew[i];
    }
    //printf("Flags: %d, %d\n", rflag, bflag);
  }

  /* output for graphical representation */
  /* Instead of using gather (which may lead to excessive memory requirements
  on the master process) each process will write its own data portion. This
  introduces a sequentialization: the hard disk can only write (efficiently)
  sequentially. Therefore, we use the following strategy:
  1. The master process writes its portion. (file creation)
  2. The master sends a signal to process 1 to start writing.
  3. Process p waites for the signal from process p-1 to arrive.
  4. Process p writes its portion to disk. (append to file)
  5. process p sends the signal to process p+1 (if it exists).
  */

  const char *filename = "ufile.txt";
  double a = 1.0;
  FILE *f;
  if (p == 0) {
    f = fopen(filename, "w");

    for (size_t i =0; i < I; i++) {
        fprintf(f, "%f ", unew[i]);
    }
    fclose(f);

    MPI_Send(&a, 1, MPI_DOUBLE, p+1, tag, MPI_COMM_WORLD);

  }else{

    MPI_Recv(&a, 1, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD, &status);

    f = fopen(filename, "a");

    for (size_t i =0; i < I; i++) {
        fprintf(f, "%f ", unew[i]);
    }
    fclose(f);
    if (p != P-1) {
      MPI_Send(&a, 1, MPI_DOUBLE, p+1, tag, MPI_COMM_WORLD);
    }
  }

  /* That's it */
  MPI_Finalize();
  exit(0);
}