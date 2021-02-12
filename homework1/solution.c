#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <complex.h>

#define WIDTH 2048
#define HEIGHT 2048

/**
 * Calculate value of each pixel - function from lecture slide
 */
int cal_pixel(complex c, double b, int N)
{
    complex z = 0 + I * 0;
    int count = 0;
    while ((cabs(z) < b) && (count < N))
    {
        z = z * z + c;
        count++;
    }
    return count;
}

int main(int argc, char **argv)
{
    int rank, nprocs;
    double b = 2.0;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Max iteration, corresponding to colors from 0 - 255
    unsigned int N = 255; 

    // Check if the number of columns could be evenly divisible by
    // the number of processor
    if (WIDTH % nprocs != 0)
    {
        fprintf(stderr, "# columns not divisible by the number of processors\n");
        MPI_Finalize();
        exit(1);
    }

    // Width of the portion of data for each processor
    int col_width = WIDTH / nprocs;
    unsigned char column[WIDTH * col_width];

    // Major work of each slave
    double zoom = 2.0; // 0.05
    double dx = zoom * b / (WIDTH - 1);
    double dy = zoom * b / (HEIGHT - 1);

    double x_offset = -1.0; // 1.302;
    double y_offset = -1.0; // 1.542;

    double complex c;
    double dreal, dimag;
    int n = 0;
    for (int x = rank * col_width; x < col_width * (rank + 1); x++)
    {
        dreal = x * dx - b + x_offset;
        for (int y = 0; y < HEIGHT; y++)
        {
            dimag = y * dy - b + y_offset;
            c = dreal + I * dimag;
            partition[n] = cal_pixel(c, b, N);
            ++n;
        }
    }

    printf("Processor %d finished.\n", rank);

    unsigned char data[HEIGHT * WIDTH];
    MPI_Gather(partition, col_width * HEIGHT, MPI_UNSIGNED_CHAR, &data,
        col_width * HEIGHT, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

    if (rank == 0) // Master
    {
        FILE *fp;
        fp = fopen(argv[1], "w");
        for (int x = 0; x < WIDTH; x++)
        {
            for (int y = 0; y < HEIGHT; y++)
            {
                fprintf(fp, "%hhu ", data[x * WIDTH + y]);
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
    }
    MPI_Finalize();
    return 0;
}