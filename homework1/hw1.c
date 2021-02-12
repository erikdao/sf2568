/**
 * Mandelbrot set parallel computation using MPI
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <mpi.h>

#define PIXELS 2048 

/* Structure definition for complex numbers */
// typedef struct {
//     double real, imag;
// } complex;

int master_code(int nprocs, int width, int height,
    double real_min, double real_max,
    double imag_min, double imag_max,
    int maxiter);

int slave_code(int myID, int width, int height, double real_min,
    double real_max, double imag_min, double imag_max, int maxiter
);

int cal_pixel(complex c, double b) {
    int max_iter = 256;
    int count = 1;
    complex z = 0 + I*0;

    while ((cabs(z) < b) && (count < max_iter)) {
        z = z*z + c;
        count++;
    }
    return count;
}

int main(int argc, char **argv)
{
    int nprocs;
    int myid;
    int maxiter;

    double real_min = -2;
    double real_max = 2;
    double imag_min = -2;
    double imag_max = 2;

    // Dimensions of display window
    int width = PIXELS;
    int height = PIXELS;
    /* Initialize for MPI */
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
        fprintf(stderr, "MPI initialization error\n");
        exit(EXIT_FAILURE);
    }
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if (nprocs < 2) {
        fprintf(stderr, "Number of processes must be at least 2\n");
        MPI_Finalize(); exit(EXIT_FAILURE);
    }

    if (myid == 0) {
        master_code(nprocs-1, width, height, real_min,
            real_max, imag_min, imag_max, maxiter);
    } else {
        slave_code(myid, width, height, real_min,
            real_max, imag_min, imag_max, maxiter);
    }

    MPI_Finalize();
    return 0;
}

int master_code(int nprocs, int width, int height,
    double real_min, double real_max,
    double imag_min, double imag_max,
    int maxiter)
{
    long start_msg[2];
    int first_row, this_row;
    int rows_per_worker;
    MPI_Status status;
    FILE *fp;
    long *data_msg = malloc((1+width) * sizeof(*data_msg));
    // Send each worker a "start" message with first row to work on
    // and number of rows

    // Calculate which rows each worker should work on
    first_row = 0;
    rows_per_worker = (int)height / nprocs;
    // long *data_msg = malloc((rows_per_worker * width) * sizeof(long));

    for (int p=0; p < nprocs; ++p) {
        start_msg[0] = first_row;
        start_msg[1] = rows_per_worker;
        MPI_Send(start_msg, 2, MPI_LONG, p+1, 2, MPI_COMM_WORLD);
        first_row += rows_per_worker;
    }

    
    // Receive the results from workers
    fp = fopen("mandelbrot_colors.txt", "w");
    for (int row = 0; row < height; ++row) {
        MPI_Recv(data_msg, width+1, MPI_LONG, MPI_ANY_SOURCE,
            100, MPI_COMM_WORLD, &status);
        for (int col = 0; col < width; ++col) {
            fprintf(fp, "%hhu ", (unsigned char)data_msg[col+1]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);

    free(data_msg);
    MPI_Finalize();
    return 0;
}

int slave_code(int myID, int width, int height, double real_min,
    double real_max, double imag_min, double imag_max, int maxiter
)
{
    long start_msg[2];
    MPI_Status status;
    int first_row, rows;
    long min_color = 0, max_color = 255;
    double scale_real, scale_imag, scale_color;

    // Receive "start" message and extract data
    MPI_Recv(start_msg, 2, MPI_LONG, 0, 2, MPI_COMM_WORLD, &status);
    first_row = start_msg[0];
    rows = start_msg[1];
    long *data_msg = malloc((1+width) * sizeof(*data_msg));

    // Compute factors to scale computational region to window
    scale_real = (double) (real_max - real_min) / (double) width;
    scale_imag = (double) (imag_max - imag_min) / (double) height;

    // Compute facot for color scaling
    scale_color = (double) (max_color - min_color) / (double) (maxiter - 1);

    // Calculate points and send back to master
    for (int row = first_row; row < first_row + rows; ++row) {
        data_msg[0] = 0;
        for (int col = 0; col < width; ++col) {
            // complex z, c;
            // z.real = z.imag = 0;

            // // Scale display coordinates to actual region
            // c.real = real_min + ((double) col * scale_real);
            // c.imag = imag_min + ((double) (height-1-row) * scale_imag);

            // int k = 0;
            // double lengthsq, temp;
            // do {
            //     temp = z.real*z.real - z.imag*z.imag + c.real;
            //     z.imag = 2*z.real*z.imag + c.imag;
            //     z.real = temp;
            //     lengthsq = z.real*z.real + z.imag*z.imag;
            //     ++k;
            // } while (lengthsq < (2*2) && k < maxiter);
            double dreal = real_min + ((double) col * scale_real);
            double dimag = imag_min + ((double) (height-1-row) * scale_imag);
            complex c = dreal + I*dimag;
            int k = cal_pixel(c, 2.0);
            // Scale color and store
            long color = (long)((k-1) * scale_color) + min_color;
            data_msg[col+1] = color;
        }
        MPI_Send(data_msg, 1 + width, MPI_LONG, 0, 100, MPI_COMM_WORLD);
    }
    printf("Node %d finished\n", myID);
    free(data_msg);
    return 0;
}