#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

int gcd(int a, int b) {
    while(a != b) {
        if(a > b) {
            a -= b;
        }
        else {
            b -= a;
        }
    }
    return a;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if (world_size < 2) {
        fprintf(stderr, "World size must be greater than 1\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int dims[1] = {world_size};
    int periods[1] = {1};
    int reorder = 0;

    MPI_Comm list_communicator;
    MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, reorder, &list_communicator);
    int my_rank;
    MPI_Comm_rank(list_communicator, &my_rank);

    int my_number = atoi(argv[my_rank + 1]);
    int got_number;

    int steps = log2(world_size);
    int left;
    int right;
    MPI_Request request;

    for(int i = 0; i < steps; i++) {
        MPI_Cart_shift(list_communicator, 0, i+1, &left, &right);

        MPI_Isend(&my_number, 1, MPI_INT, right, 0, list_communicator, &request);
        MPI_Recv(&got_number, 1, MPI_INT, left, 0, list_communicator, MPI_STATUS_IGNORE);
        
        printf("Processes: %d <- %d, Numbers: %d | %d, Step: %d\n", my_rank, left, my_number, got_number, i+1);
        my_number = gcd(my_number, got_number);
        printf("Result %d : %d\n", my_rank, my_number);
    }

    printf("Greatest common divisor from all the numbers:%d\n", my_number);
    MPI_Comm_free(&list_communicator);
    MPI_Finalize();
    return EXIT_SUCCESS;
}