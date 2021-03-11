#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int compare(const void *, const void *);

int main(int argc, char *argv[]) {
    if (argc < 2) {
        printf("N needs to be given\n");
	exit(1);
    }

    // Local variables
    int N = atoi(argv[1]);
    int I = N;
    // Data is load-balanced linearly distributed into processes
    // int I = (N + size - rank - 1) / size;
    srandom(1);
    double *x = malloc(sizeof(double) * I);

    for (int i = 0; i < I; i++) {
        x[i] = ((double) random()) / RAND_MAX;
    }

    clock_t begin = clock();
    // Local sort
    qsort(x, I, sizeof(double), compare);
    clock_t end = clock();

    printf("Time spent %f\n", (double)(end - begin)/CLOCKS_PER_SEC);

    return 0;
}

int compare(const void *x, const void *y) {
    double xx = *(double*)x, yy = *(double*)y;
    if (xx < yy) return -1;
    if (xx > yy) return 1;
    return 0;
}

