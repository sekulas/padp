#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define BUFFOR_SIZE 80

double vec_sum(double* vec, int n) {
    double sum = 0.0f;
    for(int i = 0; i < n; i++) {
        sum += vec[i];
    }
    return sum;
}

int main(int argc, char** argv) {
    MPI_Init(NULL, NULL);

    int word_size;
    MPI_Comm_size(MPI_COMM_WORLD, &word_size);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int i;
    int n;
    double *vec;
    double *res;
    double sum;

    if(rank == 0) {
        FILE *f = fopen("vector.dat", "r");
        char buffor[BUFFOR_SIZE+1];

        fgets(buffor, BUFFOR_SIZE, f);
        n = atoi(buffor);

        if(n % word_size != 0) {
            fprintf(stderr, "Amount of vec elems mod(processes) should be 0. processes: %d, vec elems: %d\n", word_size, n);
            fclose(f);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        vec = malloc(sizeof(double) * n); 
        res = malloc(sizeof(double) * word_size); 

        printf("Vector has %d elems.\n", n);

        for(i = 0; i < n; i++) {
            fgets(buffor, BUFFOR_SIZE, f);
            vec[i] = atof(buffor);
        }

        fclose(f);

        printf("Vector sum before parallel calculations: %f\n", vec_sum(vec, n));
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int elems_per_proc = n / word_size;
    double *got_vec = malloc(sizeof(double) * elems_per_proc);
    double partial_sum = 0.0f;

    MPI_Scatter(vec, elems_per_proc, MPI_DOUBLE, 
                got_vec, elems_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    partial_sum += vec_sum(got_vec, elems_per_proc);

    MPI_Gather(&partial_sum, 1, MPI_DOUBLE, 
                res, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if(rank == 0) {
        printf("Parallel sum:%f", vec_sum(res, word_size));
        free(vec);
        free(res);
    }

    free(got_vec);
    MPI_Finalize();
    return EXIT_SUCCESS;
}