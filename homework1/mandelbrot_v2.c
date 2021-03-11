#include "stdio.h"
#include "mpi.h"

#define W 2048 // Screen width
#define H 2048 // Screen height

int main(int argc, char const *argv[])
{
    int rank, num_proc, tag, rc;

    rc = MPI_Init(&argc, &argv);
    rc = MPI_Comm_Size(MPI_COMM_WORLD, &num_proc); // Get number of processors
    rc = MPI_Comm_Rank(MPI_COMM_WORLD, &rank);
    return 0;
}
