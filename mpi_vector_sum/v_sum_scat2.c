#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define BUFFOR_SIZE 80

int main(int argc, char** argv) {
    MPI_Init(NULL, NULL);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int sent_vec_n;
    int elems_per_vec;
    int* main_arr;
    int* got_arr;
    int sum = 0;
    int global_sum;

    if (rank == 0) {
        srand(time(0));
        int n = (double) rand() / RAND_MAX * (20 - 5) + 5;

        int rest = n % world_size;
        sent_vec_n = n + (rest > 0 ? (world_size - rest) : 0);

        printf("N: %d\nSent n size: %d\n", n, sent_vec_n);

        main_arr = (int*) malloc(sizeof(int) * sent_vec_n);

        for(int i = 0; i < n; i++) {
            main_arr[i] = 1;
        }
        for(int i = n; i < sent_vec_n; i++) {
            main_arr[i] = 0;
        }
    }

    MPI_Bcast(&sent_vec_n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    elems_per_vec = sent_vec_n / world_size;

    got_arr = (int*) malloc(sizeof(int) * elems_per_vec);

    MPI_Scatter(main_arr, elems_per_vec, MPI_INT, got_arr, elems_per_vec, MPI_INT, 0, MPI_COMM_WORLD);

    for(int i = 0; i < elems_per_vec; i++) {
        sum += got_arr[i];
    }

    MPI_Reduce(&sum, &global_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("GLOBAL SUM: %d", global_sum);
        free(main_arr);
    }

    free(got_arr);
    MPI_Finalize();
}