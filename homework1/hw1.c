/**
 * Mandelbrot set parallel computation using MPI
*/
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define PIXELS 2048

// Structure representing complex number
typedef struct complex
{
    float real;
    float imag;
} Complex;

int cal_pixel(Complex c, int b) {
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
    } while ((lengthsq < (b * b)) && (count < max_iter));
    return count;
}

int main(int argc, char **argv)
{
    int nprocs;
    int myid;

    // Dimensions of display window
    int width = PIXELS;
    int height = PIXELS;
    int *data = malloc((width+1) * sizeof(*data));
    
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if (nprocs < 2) {
        fprintf(stderr, "Number of processes must be at least 2\n");
        MPI_Finalize();
        exit(-1);
    }

    int first_row, number_of_rows;

    if (myid == 0) { // Master
        // Send each workder a message with first row to work on and number of rows
        first_row = 0;
        number_of_rows = height / nprocs;
        
    }
    MPI_Finalize();
    return 0;
}
