/**
 * Mandelbrot set parallel computation using MPI
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <mpi.h>

#define PIXELS 2048

// Structure representing complex number
// typedef struct complex
// {
//     float real;
//     float imag;
// } Complex;

// int cal_pixel(Complex c, int b) {
//     int count, max_iter;
//     Complex z;
//     float temp, lengthsq, abs_z;

//     // Init
//     max_iter = 256;
//     z.real = 0;
//     z.imag = 0;
//     count = 0;

//     do {
//         temp = z.real * z.real - z.imag * z.imag + c.real;
//         z.imag = 2 * z.real * z.imag + c.imag;
//         z.real = temp;
//         // lengthsq = z.real * z.real + z.imag * z.imag;
//         abs_z = pow((z.real * z.real + z.imag * z.imag), 0.5);
//         // printf("abs_z: %f\n", abs_z);
//         count++;
//     } while ((abs_z < b) && (count < max_iter));
//     return count;
// }

int pixel_value(double complex c, float b, int N) {
    int count = 1;
    double complex z = 0 + I*0;

    while ((cabs(z) <= b) && (count < N)) {
        z = z*z + c;
        count++;
    }
    return count;
}

int main(int argc, char **argv)
{
    int nprocs;
    int rank;
    int rc;

    // Dimensions of display window
    int width = PIXELS;
    int height = PIXELS;
    int *data_msg = malloc((width+1) * sizeof(*data_msg));
    FILE *fp;

    MPI_Status status;

    rc = MPI_Init(&argc, &argv);
    rc = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    rc = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    printf("Rank %d", rank);

    int first_row, number_of_rows;
    int start_msg[2];
    int k; // process iterator
    int row, col;
    int b = 2;
    double dx = (double) 2 * b / (width-1);
    double dy = (double) 2 * b / (height-1);
    double dreal, dimag;
    int this_row;

    if (rank == 0) { // Master
        // Send each slave a message with first row to work on and number of rows
        first_row = 0;
        number_of_rows = height / nprocs;
        printf("first_row: %d, number_of_rows: %d", first_row, number_of_rows);
        for (k = 1; k < nprocs; k++) {
            start_msg[0] = first_row;
            start_msg[1] = number_of_rows;
            MPI_Send(start_msg, 2, MPI_INTEGER, k, 2, MPI_COMM_WORLD);
            first_row += number_of_rows;
        }

        // Receive results from slaves
        fp = fopen("output.txt", "w");
        for (row = 0; row < height; ++row) {
            MPI_Recv(data_msg, width+1, MPI_INTEGER, MPI_ANY_SOURCE, 100, MPI_COMM_WORLD, &status);
            this_row = data_msg[0];
            printf("this_row: %d", this_row);
            for (col = 0; col < width; ++col) {
                fprintf(fp, "%hhu ", data_msg[col+1]);
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
    } else { // Slave
        // Receive "start" message and extract data
        MPI_Recv(start_msg, 2, MPI_INTEGER, 0, 2, MPI_COMM_WORLD, &status);
        first_row = start_msg[0];
        number_of_rows = start_msg[1];

        // Calculate points and send back to master a row at a time
        data_msg = malloc((width+1) * sizeof(*data_msg));
        for (row = first_row; row < first_row + number_of_rows; ++row) {
            data_msg[0] = row;
            dreal = row * dx - b;
            for (col = 0; col < width; col++) {
                dimag = col * dy - b;
                double complex c = dreal + I*dimag;
                data_msg[col+1] = pixel_value(c, b, 256);
            }
        }

        MPI_Send(data_msg, width+1, MPI_INTEGER, 0, 100, MPI_COMM_WORLD);
        printf("Node %d finalized, gathering data. \n", rank);
    }
    MPI_Finalize();
    free(data_msg);
    return 0;
}
