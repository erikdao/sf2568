#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void my_bcast(
    void* data, int count, MPI_Datatype dtype, int root,
    MPI_Comm communicator
) {
    int rank, size;
    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator, &size);

    if (rank == root) {
        // Root process, send data to everyone
        int i = 0;
        for (i = 0; i < size; i++) {
            if (i != rank) {
                MPI_Send(data, count, dtype, i, 0, communicator);
            }
        }
    } else {
        // Receiver process
        MPI_Recv(data, count, dtype, root, 0, communicator, MPI_STATUS_IGNORE);
    }
}

int main(int argc, char** argv) {
    MPI_Init(NULL, NULL);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int data;
    if (world_rank == 0) {
        data = 100;
        printf("Process 0 broadcasting data %d\n", data);
        my_bcast(&data, 1, MPI_INT, 0, MPI_COMM_WORLD);
    } else {
        my_bcast(&data, 1, MPI_INT, 0, MPI_COMM_WORLD);
        printf("Process %d received data %d from root\n", world_rank, data);
    }
    MPI_Finalize();
}