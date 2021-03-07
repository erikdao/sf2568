#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

/**
 * Returns -1 if a < b
 *          1 if a > b
 *          0 if a = b
 */
int cmpfunc(const void * x, const void * y) {
    double xx = *(double*)x, yy = *(double*)y;
    if (xx < yy) return -1;
    if (xx > yy) return 1;
    return 0;
}


int main(int argc, char **argv) {
    int P, myrank, N, I;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    /* Find problem size N from command line */
    if (argc < 2) {
        printf("No size N given");
        exit(1);
    }
    N = atoi(argv[1]);

    // Local size
    I = (N + P - myrank - 1) / P; 
    srandom(myrank + 1);

    // Generate data
    double *x = malloc(sizeof(double) * I);
    for (int i = 0; i < I; i++) {
        x[i] = ((double) random()) / RAND_MAX;
    }
    // Local sorting
    qsort(x, I, sizeof(double), cmpfunc);

    // Global sorting phase
    for (int step = 0; step < P / 2; step++) {

    }
    MPI_Finalize();
    return 0;
}