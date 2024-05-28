#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

const int AMOUNT_OF_DATA = 1000000;
const int PACKET_SIZE = 1000;
const int ITERATIONS = AMOUNT_OF_DATA / PACKET_SIZE;

int main() {
    MPI_Init(NULL, NULL);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if(world_size < 2) {
        fprintf(stderr, "World size must be greater than 1\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int dims[1] = {world_size};
    int periods[1] = {0};
    int reorder = 0;

    MPI_Comm pipe;
    MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, reorder, &pipe);

    int my_rank;
    MPI_Comm_rank(pipe, &my_rank);

    int got_packet[PACKET_SIZE];

    int left, right;
    MPI_Cart_shift(pipe, 0, 1, &left, &right);

    if(my_rank == 0) {
        int data[AMOUNT_OF_DATA];

        for (int i = 0; i < AMOUNT_OF_DATA; i++) {
            data[i] = 1;
        }

        for(int i = 0; i < ITERATIONS; i++) {
            MPI_Send(&data[i*PACKET_SIZE], PACKET_SIZE, MPI_INT, right, 0, pipe);
        }
    }
    else {
        if(right != MPI_PROC_NULL){
            for(int i = 0; i < ITERATIONS; i++) {
                MPI_Recv(&got_packet, PACKET_SIZE, MPI_INT, left, 0, pipe, MPI_STATUS_IGNORE);
                MPI_Send(&got_packet, PACKET_SIZE, MPI_INT, right, 0, pipe);

            }
        }
        else {
            int sum = 0;

            for(int i = 0; i < ITERATIONS; i++) {
                MPI_Recv(&got_packet, PACKET_SIZE, MPI_INT, left, 0, pipe, MPI_STATUS_IGNORE);

                for(int j = 0; j < PACKET_SIZE; j++){
                    sum += got_packet[j];
                }
            }

            printf("Final sum: %d", sum);
        }
    }

    MPI_Comm_free(&pipe);
    MPI_Finalize();
    return EXIT_SUCCESS;
}