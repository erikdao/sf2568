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

void printArray(double *x, int len){
    for( int i = 0; i<len; i++){
        printf("%f ", x[i]);
    }
    printf("\n");
}

void mergeArrays(double *arr1, double *arr2, int n1, int n2, double *arr3)
{
    int i = 0, j = 0, k = 0;

    // Traverse both array
    while (i<n1 && j <n2)
    {
        // Check if current element of first
        // array is smaller than current element
        // of second array. If yes, store first
        // array element and increment first array
        // index. Otherwise do same with second array
        if (arr1[i] < arr2[j])
            arr3[k++] = arr1[i++];
        else
            arr3[k++] = arr2[j++];
    }

    // Store remaining elements of first array
    while (i < n1)
        arr3[k++] = arr1[i++];

    // Store remaining elements of second array
    while (j < n2)
        arr3[k++] = arr2[j++];
    printf("k = %d\n", k);
}

int main(int argc, char **argv) {
    int P, myrank, N, I;

    MPI_Init(&argc, &argv);
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    /* Find problem size N from command line */
    if (argc < 2) {
        printf("No size N given");
        exit(1);
    }
    N = atoi(argv[1]);

    // Local size
    //I = (N + P - myrank - 1) / P;
    I = N/P;
    if(N%P != 0)
        {I += 1;}
    double *x = malloc(sizeof(double) * I);
    
    // Generate data
    //fix from here
    /*
    if(N%P != 0){
        I += 1;
        srandom(myrank +1);
        double *x = malloc(sizeof(double) * (I+1));
        for (int i = 0; i < I; i++) {
            x[i] = ((double) random()) / RAND_MAX;
        }
        if (myrank >= N%P){
            x[I] = -1.0;
        }
    }
    else{
        I = N/P;
        srandom(myrank + 1);
        double *x = malloc(sizeof(double) * I);
        for (int i = 0; i < I; i++) {
            x[i] = ((double) random()) / RAND_MAX;
        }
    }
    */
    double *a = malloc(sizeof(double) * I);
    // double *result = malloc(length_x * sizeof(double) * 2);
    srandom(myrank +1);
    for (int i = 0; i < I; i++) {
        x[i] = ((double) random()) / RAND_MAX;
        }

    if(N%P != 0){
        printf("modulo: %d\n",N%P);
        printf("myrank: %d\n",myrank);
        if (myrank >= N%P){
            x[I-1] = -1.0;
        }
    }

    // Local sorting
    qsort(x, I, sizeof(double), cmpfunc);

    // Global sorting phase
    int evenproc = ((myrank % 2) == 0);
    int evenphase = 1;

    // Initialize the result array, its size will be changed later
    double *result = malloc(I * sizeof(double));
    int size_before,size_after; 
    for (int step = 0; step < P/2; step++) {
	printf("Rank: %d; step: %d; evenphase: %d\n", myrank, step, evenphase);
        size_before = I * pow(2, step); 
        size_after = I * pow(2, step + 1);
	printf("size_before: %d; size_after: %d\n", size_before, size_after);
	printf("x\t"); printArray(x, size_before);

        if (evenproc && evenphase){
            MPI_Recv(a,size_before,MPI_DOUBLE,myrank+1,100,MPI_COMM_WORLD,&status);
            MPI_Send(x,size_before,MPI_DOUBLE,myrank+1,100,MPI_COMM_WORLD);
        }
        else if(evenproc && !evenphase){
            if(myrank >= 2){
                MPI_Recv(a,size_before,MPI_DOUBLE,myrank-1,100,MPI_COMM_WORLD,&status);
                MPI_Send(x,size_before,MPI_DOUBLE,myrank-1,100,MPI_COMM_WORLD);
            }
        }
        else if(!evenproc && evenphase){
            MPI_Send(x,size_before,MPI_DOUBLE,myrank-1,100,MPI_COMM_WORLD);
            MPI_Recv(a,size_before,MPI_DOUBLE,myrank-1,100,MPI_COMM_WORLD,&status);
        }
        else if(!evenproc && !evenphase){
            if(myrank <= P - 3){
                MPI_Send(x,size_before,MPI_DOUBLE,myrank+1,100,MPI_COMM_WORLD);
                MPI_Recv(a,size_before,MPI_DOUBLE,myrank+1,100,MPI_COMM_WORLD,&status);
            }
        }
	result = realloc(result, size_after * sizeof(double));
    mergeArrays(x,a,size_before,size_before,result);
	printf("Merged: "); printArray(result, size_after);
        x = realloc(x, sizeof(double) * size_after);
        a = realloc(a, sizeof(double) * size_after);
        x = result;
        evenphase = !evenphase;
	printf("-------------------------------\n");
    }
    int size_final = size_after - N%P;
    double *final = malloc(sizeof(double) * (size_final) );
    for (int i = 0; i < size_final; i++) {
        final[i] = result[i + (N%P)];
    }
    printf("Final: "); printArray(final, size_final);
    // free(x);
    // free(a);
    // free(result);
    MPI_Finalize();
    return 0;
}

