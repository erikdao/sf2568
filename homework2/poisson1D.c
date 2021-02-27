/* Reaction-diffusion equation in 1D
 * Solution by Jacobi iteration
 * simple MPI implementation
 *
 * C Michael Hanke 2006-12-12
 */
#include <stdlib.h>
#include <stdio.h>

#define MIN(a,b) ((a) < (b) ? (a) : (b))

/* Use MPI */
#include "mpi.h"

/* define problem to be solved */
#define N 1000   /* number of inner grid points */
#define SMX 1000000 /* number of iterations */

/* implement coefficient functions */
extern double r(const double x);
extern double f(const double x); 

double r(const double x) {
    // Change here before submission: -x
    // return x - 1;
    // return -x;
    return x - 1;
}

double f(const double x) {
    // return -2 + x * x * (x - 1);
    // return 2 + (x - 1) * x * (x - 1);
    return 6*x - 3 + (x-1)*x*(x-0.5)*(x-1);
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
    int P, p, tag, L, R, I, n, color, step, i;
    double h, *u, *unew, *ff, *rr, x_n;
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
    color = p % 2;  // 0: red, 1: black
    printf("Node %d, color: %d\n", p, color);
    /* Compute local indices for data distribution */
    L = N / P;
    R = N % P;
    I = (N + P - p - 1) / P;  // Number of local elements
    h = 1.0 / (N + 1);

    printf("Processor %d, L = %d; R = %d; I = %d; h = %f\n", p, L, R, I, h);
    /* arrays */
    unew = (double *) malloc(I*sizeof(double));
    /* Note: The following allocation includes additionally:
    - boundary conditins are set to zero
    - the initial guess is set to zero */
    u = (double *) calloc(I+2, sizeof(double));

    ff = (double *) malloc(I * sizeof(double));
    rr = (double *) malloc(I * sizeof(double));

    for (i = 0; i < I; i++) {
        n = p * L + MIN(p, R) + i;  // Global index for given (p, i)
        x_n = n * h;
        ff[i] = f(x_n);
        rr[i] = r(x_n);
    }

    /* Jacobi iteration */
    for (step = 0; step < SMX; step++) {
        if (step % 100000 == 0) {
            printf("Jacobi, processor %d, iter %d\n", p, step);
        }
        /* RB communication of overlap */
        if (color == 0) { // red
            MPI_Send(&u[I], 1, MPI_DOUBLE, p+1, tag, MPI_COMM_WORLD);
            MPI_Recv(&u[I+1], 1, MPI_DOUBLE, p+1, tag, MPI_COMM_WORLD, &status);
            if (p != 0) {
                MPI_Send(&u[1], 1, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD);
                MPI_Recv(&u[0], 1, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD, &status);
            }
        } else { // black
            MPI_Recv(&u[0], 1, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD, &status);
            MPI_Send(&u[1], 1, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD);
            if (p != P-1) {
                MPI_Recv(&u[I+1], 1, MPI_DOUBLE, p+1, tag, MPI_COMM_WORLD, &status);
                MPI_Send(&u[I], 1, MPI_DOUBLE, p+1, tag, MPI_COMM_WORLD);
            }
        }
        /* local iteration step */
        for (i = 0; i < I; i++) {
            unew[i] = (u[i]+u[i+2]-h*h*ff[i])/(2.0-h*h*rr[i]);
            // unew[i] = (u[i-1] + u[i+1] - h*h*ff[i]) / (2.0 - h*h*rr[i]);
        }
        for (i = 0; i < I; i++) {
            u[i+1] = unew[i];
        }
    }

    /* output for graphical representation */
    /* Instead of using gather (which may lead to excessive memory requirements
    on the master process) each process will write its own data portion. This
    introduces a sequentialization: the hard disk can only write (efficiently)
    sequentially. Therefore, we use the following strategy: */
    FILE *f;
    double signal;
    // 1. The master process writes its portion. (file creation)
    if (p == 0) {
        f = fopen("possion.txt", "w");
        for (i = 0; i < I; i++) {
            fprintf(f, "%f ", unew[i]);
        }
        fclose(f);
        // 2. The master sends a signal to process 1 to start writing.
        signal = 1.0;
        MPI_Send(&signal, 1, MPI_DOUBLE, 1, tag, MPI_COMM_WORLD);
        printf("Master processor wrote and sent signal\n");
    } else {
        // 3. Process p waites for the signal from process p-1 to arrive.
        MPI_Recv(&signal, 1, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD, &status);
        printf("Received signal %f from processor %d\n", signal, p-1);
        if (signal == 1.) {
            // 4. Process p writes its portion to disk. (append to file)
            f = fopen("possion.txt", "a");
            for (i = 0; i < I; i++) {
                fprintf(f, "%f ", unew[i]);
            }
            fclose(f);
            printf("Processor %d wrote its portion\n", p);
            // 5. process p sends the signal to process p+1 (if it exists).
            if (p != P - 1) {
                MPI_Send(&signal, 1, MPI_DOUBLE, p+1, tag, MPI_COMM_WORLD);
                printf("Processor %d sent signal to processor %d\n", p, p+1);
            }
        }
    }

    /* That's it */
    MPI_Finalize();
    exit(0);
}

