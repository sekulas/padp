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

    integration_data_t integration_data;
    double partial_integral;

    MPI_Status status;
    MPI_Request request;

    switch(world_rank) {
        case 0:
            int segments = num_of_points - 1;
            int segments_per_process = segments / world_size;
            int rest_of_segments = segments % world_size;
            double dx = (end - begin) / segments;
            int segments_for_process;

            for(int i = 0; i < world_size - 1; i++) {
                segments_for_process = (segments_per_process + (rest_of_segments > 0 ? 1 : 0));
                set_integration_data(&integration_data, &begin, dx, segments_for_process);
                rest_of_segments--;
                printf("SENDING: %d, Begin: %lf, End: %lf\n", i, integration_data.begin, integration_data.end);
                MPI_Isend(&integration_data, 1, mpi_integration_data_t, i + 1, 1, MPI_COMM_WORLD, &request);
            }

            segments_for_process = (segments_per_process + (rest_of_segments > 0 ? 1 : 0));
            set_integration_data(&integration_data, &begin, dx, segments_for_process);
            double global_integral = integrate(x2, integration_data.begin, 
                                                integration_data.end, integration_data.segments);
            printf("MOTHER: %d, Begin: %lf, End: %lf\n", world_rank, integration_data.begin, integration_data.end);

            MPI_Status* statuses = malloc((world_size - 1) * sizeof(MPI_Status));
            MPI_Request* requests = malloc((world_size - 1) * sizeof(MPI_Request));
            double* partial_integrals = malloc((world_size - 1) * sizeof(double));

            for(int i = 1; i < world_size; i++) {
                MPI_Irecv(&partial_integrals[i-1], 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &requests[i-1]);
            }

            MPI_Waitall(world_size-1, requests, statuses);

            for(int i = 0; i < world_size - 1; i++) {
                global_integral += partial_integrals[i];
            }

            printf("Integral: %lf\n", global_integral);
            free(statuses);
            free(requests);
            free(partial_integrals);
            break;
        default:
            MPI_Irecv(&integration_data, 1, mpi_integration_data_t, 0, 1, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, &status);
            printf("My rank: %d, Begin: %lf, End: %lf\n", world_rank, integration_data.begin, integration_data.end);
            partial_integral = integrate(x2, integration_data.begin, 
                                            integration_data.end, integration_data.segments);
            
            MPI_Isend(&partial_integral, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &request);
            break;
    }
    
    MPI_Type_free(&mpi_integration_data_t);
    MPI_Finalize();
}