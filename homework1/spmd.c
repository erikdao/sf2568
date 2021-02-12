#include "stdio.h"
#include "stdlib.h"
#include "complex.h"
#include "mpi.h"


int cal_pixel(complex c, double b, int N) {
    complex z = 0 + I*0;
    int count = 1;
    while ((cabs(z) < b) && (count < N)) {
        z = z*z + c;
        count++;
    }
    return count;
}

int main(int argc, char **argv)
{
    int rank, nprocs, rc;
    int width = 2048;
    int height = 2048;
    int b = 2.0;  // Bound
    int N = 255;  // Max iteration - corresponding to 0 - 255 colors
    FILE *fp;
    MPI_Status status;

    rc = MPI_Init(&argc, &argv);
    rc = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    rc = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (height % nprocs != 0) {
        fprintf(stderr, "Height %d is not divisible by number of processes %d\n", height, nprocs);
        MPI_Finalize();
        exit(1);
    }
    int partition_width = width / nprocs;
    unsigned char partition[height * partition_width];

    if (rank == 0) {
        unsigned char data[width*height];

        for (int k = 1; k < nprocs; k++) {
            MPI_Recv(partition, height * partition_width, MPI_UNSIGNED_CHAR,
                MPI_ANY_SOURCE, 100, MPI_COMM_WORLD, &status);
            printf("Received message (data) from process %d\n", k);
            for (int i = 0; i < height * partition_width; i++) {
                if (k == 1) {
                    data[i] = partition[i];
                } else {
                    data[(k-1)*height*partition_width + i] = partition[i];
                }
            }
        }

        fp = fopen("mandelbrot.txt", "w");
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                fprintf(fp, "%hhu ", data[x * width + y]);
            }
            fprintf(fp, "\n");
        }
        fclose(fp);

    } else { // Slave
        double dreal, dimag;
        double dx = (double) (2 * b / (width - 1));
        double dy = (double) (2 * b / (height - 1));
        int wp = width / nprocs;
        int hp = height;
        int xoff = rank * wp;
        int yoff = 0;

        int i = 0;
        // for (int x = 0; i < wp; x++) {
        for (int x = rank * partition_width; x < partition_width * (rank+1); x++) {
            // dreal = (x + xoff) * dx - b;
            dreal = x * dx - b;
            for (int y = 0; i < height; y++) {
                // dimag = (y + yoff) * dy - b;
                dimag = y * dy - b;
                complex c = dreal + I*dimag;
                partition[i] = cal_pixel(c, b, N);
                i++;
            }
        }
        MPI_Send(partition, height * partition_width, MPI_UNSIGNED_CHAR, 0, 100, MPI_COMM_WORLD);
        printf("Process %d finished, sent back data\n", rank);
    }
    MPI_Finalize();
    return 0;
}
