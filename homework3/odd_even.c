#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <time.h>

#define EMPTY_PLACEHOLDER -1.0

void print_array(double *arr, const int len);
int compute_neighbor(int phase, int rank, int size);
int compare(const void *, const void *);
void merge_arrays(double *, int, double *, int, double *, unsigned int);

int main(int argc, char *argv[]) {
    int size, rank;
    double time_spent = 0.0;

    if (argc < 2) {
        printf("N needs to be given\n");
	exit(1);
    }

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Status status;

    // Local variables
    int N = atoi(argv[1]);
    // Data is load-balanced linearly distributed into processes
    // int I = (N + size - rank - 1) / size;
    int I = N / size;
    if (N % size != 0)
	I += 1;

    // Data generation
    srandom(rank + 1);
    double *x = malloc(sizeof(double) * I);
    double *a = malloc(sizeof(double) * I);
    double *temp = malloc(sizeof(double) * I *2);

    for (int i = 0; i < I; i++) {
        x[i] = ((double) random()) / RAND_MAX;
    }

    if ((N % size != 0) && (rank >= N % size)) {
        x[I-1] = EMPTY_PLACEHOLDER;
    }
    clock_t begin = clock();
    // Local sort
    qsort(x, I, sizeof(double), compare);

    // Odd-even transposition
    for (int phase = 0; phase < size; phase++) {
	MPI_Barrier(MPI_COMM_WORLD);
	int neighbor = compute_neighbor(phase, rank, size);

	if (neighbor >= 0 && neighbor < size) {
            MPI_Sendrecv(x, I, MPI_DOUBLE, neighbor, phase,
			 a, I, MPI_DOUBLE, neighbor, phase,
			 MPI_COMM_WORLD, &status);

	    if (rank < neighbor) {
	        merge_arrays(x, I, a, I, temp, 1);
	    } else {
	        merge_arrays(x, I, a, I, temp, 0);
	    }
	}
    }
    clock_t end=clock();
    time_spent += (double)(end-begin)/CLOCKS_PER_SEC;
    printf("time spent is %f seconds\n", time_spent);

    // Sequentially write the sorted array
    int signal = 0;
    FILE *f;
    char fname[50];
    sprintf(fname, "sorted_array_s-%d_N-%s.txt", size, argv[1]);

    if (rank == 0) { // master process
    	f = fopen(fname, "w");
	for (int i = 0; i < I; i++) {
	    if (x[i] != EMPTY_PLACEHOLDER) {
	    	fprintf(f, "%1.10f\n", x[i]);
	    }
	}
	fclose(f);
	signal = 1;
	MPI_Send(&signal, 1, MPI_INT, rank+1, 100, MPI_COMM_WORLD);
    } else {
	MPI_Recv(&signal, 1, MPI_INT, rank-1, 100, MPI_COMM_WORLD, &status);
	if (signal == 1) {
	    f = fopen(fname, "a");
	    for (int i = 0; i < I; i++) {
		if (x[i] != EMPTY_PLACEHOLDER) {
	       	    fprintf(f, "%1.10f\n", x[i]);
		}
	    }
	    fclose(f);
	    if (rank != size - 1) {
	        MPI_Send(&signal, 1, MPI_INT, rank+1, 100, MPI_COMM_WORLD);
	    }
	}
    }

    free(x);

    MPI_Finalize();
    return 0;
}

void print_array(double *arr, const int len) {
    for (int i = 0; i < len; i++) {
        printf("%f ", arr[i]);
    }
    printf("\n");
}

int compute_neighbor(int phase, int rank, int size) {
    int neighbor;

    if (phase % 2 != 0) {  // Odd phase
	if (rank % 2 != 0) { // Odd processor
	    neighbor = rank + 1;
	} else { // Even processor
	    neighbor = rank - 1;
	}
    } else { // Even phase
	if (rank % 2 != 0) { // Odd processor
	    neighbor = rank - 1;
	} else { // Even processor
	    neighbor = rank + 1;
	}
    }

    if (neighbor < 0 || neighbor >= size) {
        neighbor = -1;
    }
    return neighbor;
}

int compare(const void *x, const void *y) {
    double xx = *(double*)x, yy = *(double*)y;
    if (xx < yy) return -1;
    if (xx > yy) return 1;
    return 0;
}

/**
 * Merge two **sorted** array `src` and `rec` into the `temp` array
 * `keep_low` decides if we want to keep the lower part in the `src`
 * or the upper part in the `src` array after this merge operation
 */
void merge_arrays(
    double *src, int len_src, double *recv, int len_recv,
    double *temp, unsigned int keep_low
) {
    int i = 0, j = 0, k = 0;

    while (i < len_src && j < len_recv) {
        if (src[i] < recv[j])
	    temp[k++] = src[i++];
	else
	    temp[k++] = recv[j++];
    }

    // Store remaining elements of first array
    while (i < len_src)
	temp[k++] = src[i++];

    // Store remaining elements of second array
    while (j < len_recv)
	temp[k++] = recv[j++];

    if (keep_low == 1) {
	// Keep the lower part of the sorted `temp` array in
	// `src` for this process
	for (int i = 0; i < len_src; i++) {
	    src[i] = temp[i];
	}
    } else if (keep_low == 0) {
	// Keep the upper part of the sorted `temp` array in
	// `src` for this process
    	for (int i = len_src, j = 0; j < len_recv; i++, j++) {
	    src[j] = temp[i];
	}
    }
}
