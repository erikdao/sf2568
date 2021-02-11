/**
 * Mandelbrot set parallel computation using MPI
*/
#include <stdio.h>
#include <mpi.h>

// Structure representing complex number
typedef struct complex
{
    float real;
    float imag;
} Complex;

int cal_pixel(Complex c) {
    int count, max_iter;
    Complex z;
    float temp, lengthsq;

    // Init
    max_iter = 256;
    z.real = 0;
    z.imag = 0;
    count = 0;

    do {
        temp = z.real * z.real - z.imag * z.imag + c.real;
        z.imag = 2 * z.real * z.imag + c.imag;
        z.real = temp;
        lengthsq = z.real * z.real + z.imag * z.imag;
        count++;
    } while ((lengthsq < 4.0) && (count < max_iter));
    return count;
}

int main(int argc, char **argv)
{
    int rank, num_proc, tac, rc, i, j;
    int *data;
    Complex c;
    MPI_Status status;

    rc = MPI_Init(&argc, &argv);
    rc = MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    rc = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    printf("Number of processors: %d\n", num_proc);
    printf("Current process: %d\n", rank);

    rc = MPI_Finalize();
    return 0;
}
