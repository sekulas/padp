#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <mpi.h>

double x2(double x) {
    return x * x;
}

double integrate(double (*func)(double), double begin, double end, int segments) {
    double dx = (end - begin) / segments;
    double sum = 0;
    double step = begin;

    for(int i = 0; i < segments; i++) {
        sum += (func(step) + func(step + dx)) * dx/2;
        step += dx;
    }

    return sum;
}

typedef struct integration_data {
    double begin;
    double end;
    int segments;
} integration_data_t;

void create_mpi_integral_data_type(MPI_Datatype * mpi_integral_data_t) {
    int block_lengths[3] = {1, 1, 1};
    MPI_Aint displacements[3];
    MPI_Datatype types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_INT};

    integration_data_t temp;

    displacements[0] = offsetof(integration_data_t, begin);
    displacements[1] = offsetof(integration_data_t, end);
    displacements[2] = offsetof(integration_data_t, segments);

    MPI_Type_struct(3, block_lengths, displacements, types, mpi_integral_data_t);
    MPI_Type_commit(mpi_integral_data_t);
}

void set_integration_data(integration_data_t* data, double* begin, double dx, int segments) {
    data->begin = *begin;
    *begin += dx * segments;
    data->end = *begin;
    data->segments = segments;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    double begin = atof(argv[1]);
    double end = atof(argv[2]);
    int num_of_points = atoi(argv[3]);

    MPI_Datatype mpi_integration_data_t;
    create_mpi_integral_data_type(&mpi_integration_data_t);

    integration_data_t integration_data_s[world_size];
    integration_data_t integration_data_r;
    double partial_integrals[world_size];
    double partial_integral;

    if(world_rank == 0) {
        int segments = num_of_points - 1;
        int segments_per_process = segments / world_size;
        int rest_of_segments = segments % world_size;
        double dx = (end - begin) / segments;
        int segments_for_process;

        for(int i = 0; i < world_size; i++) {
            segments_for_process = (segments_per_process + (rest_of_segments > 0 ? 1 : 0));
            set_integration_data(&integration_data_s[i], &begin, dx, segments_for_process);
            rest_of_segments--;
        }
    }

    MPI_Scatter(&integration_data_s, 1, mpi_integration_data_t, 
        &integration_data_r, 1, mpi_integration_data_t, 0, MPI_COMM_WORLD);

    partial_integral = integrate(x2, integration_data_r.begin, 
                                    integration_data_r.end, integration_data_r.segments);

    MPI_Gather(&partial_integral, 1, MPI_DOUBLE, &partial_integrals, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    if(world_rank == 0) {
        for (int i = 1; i < world_size; i++) {
            partial_integral += partial_integrals[i];
        }
        printf("Integral: %lf\n", partial_integral);
        MPI_Type_free(&mpi_integration_data_t);
        MPI_Finalize();
    }
    
    MPI_Finalize();
    return EXIT_FAILURE;
}