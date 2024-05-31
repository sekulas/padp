#include <stdio.h>
#include <stdlib.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <time.h>

#define BUFFOR_SIZE 80
#define SUBPR_ARR_LEN 6
#define KEY_IDX 65
#define KEY_RES 66
#define KEY_DAT 67
#define PROG_NAME "vec_sum.c"

enum ShmVecType {
    INDEXES,
    RESULTS,
    DATA
};

void* initialize_shm_vec(int vec_size, enum ShmVecType vec_type)
{  
    int key;
    int size_idx;

    switch(vec_type) {
        case INDEXES:
            key = KEY_IDX;
            size_idx = sizeof(int) * vec_size * 2;
            break;
        case RESULTS:
            key = KEY_RES;
            size_idx = sizeof(double) * vec_size;
            break;
        case DATA:
            key = KEY_DAT;
            size_idx = sizeof(double) * vec_size;
            break;
        default:
            fprintf(stderr, "Bad ShmVecType enum specified.\n");
            exit(1);
    }

	key_t key_idx = ftok(PROG_NAME, key);
	int seg_idx = shmget(key_idx, size_idx, 0666 | IPC_CREAT);
	if (seg_idx < 0) {
		fprintf(stderr, "shmget error in child(%d) process!", getpid());
		exit(1);
	}

	void* shm_vec = shmat(seg_idx, 0, 0);
	if (shm_vec == (int *)-1) {
		fprintf(stderr, "shmat error in child(%d) process!", getpid());
		exit(1);
	}

	return shm_vec;
}

double *read_from_file(FILE *f, int n)
{
	char buffor[BUFFOR_SIZE + 1];
	double *vector = (double*) initialize_shm_vec(n, DATA);
	printf("Vector has %d elements\n", n);
	for (int i = 0; i < n; i++) {
		fgets(buffor, BUFFOR_SIZE, f);
		vector[i] = atof(buffor);
	}
	fclose(f);
	return vector;
}

double child_sum(double* vector, int start, int end) {
	double sum = 0.0f;
	for (int i = start; i < end; i++) {
		sum += vector[i];
	}
	return sum;
}

void fill_indexes(int n, int subpr_n, int* indexes) {
    int elements_per_subprocess = n / subpr_n;
    int remaining_elements = n % subpr_n;
    int start = 0;
    int end = 0;
    int counter = 0;

    for (int i = 0; i < subpr_n; i++) {
        start = end;
        end = start + elements_per_subprocess + (i < remaining_elements ? 1 : 0);

        indexes[counter++] = start;
        indexes[counter++] = end;
    }
}

int main(int argc, char** argv) {
    int subpr_n = -1;
    int subpr_arr[SUBPR_ARR_LEN] = {1, 2, 4, 6, 8, 16};
    int n;
    double final_sum;

    double* times = (double*) malloc(sizeof(double) * SUBPR_ARR_LEN);

    for(int i = 0; i < SUBPR_ARR_LEN; i++) {
        subpr_n = subpr_arr[i];
        double final_sum = 0;

        pid_t pid;
        pid_t* kids = (pid_t*) malloc(sizeof(pid_t) * subpr_n);

        FILE *f = fopen("vector.dat", "r");
        if(f == NULL) {
            fprintf(stderr, "Cannot open 'vector.dat'.\n");
            return EXIT_FAILURE;
        }

        char buffor[BUFFOR_SIZE + 1];
        fgets(buffor, BUFFOR_SIZE, f);
        n = atoi(buffor);

        clock_t timer;
        timer = clock();

        int* indexes = (int*) initialize_shm_vec(subpr_n, INDEXES);
        double* results = (double*) initialize_shm_vec(subpr_n, RESULTS);
        double *vector;
        double sum;

        for(int j = 0; j < subpr_n; j++) {
            switch(pid = fork()) {
                case -1:
                    fprintf(stderr, "Error in fork.\n");
                    return EXIT_FAILURE;
                case 0:
                    sigset_t mask;
                    struct sigaction onsig;
                    sigfillset(&mask);
                    sigdelset(&mask, SIGUSR1);
                    onsig.sa_handler = NULL;
                    onsig.sa_mask = mask;
                    onsig.sa_flags = SA_SIGINFO;
                    sigaction(SIGUSR1, &onsig, NULL);
                    
                    vector = (double*) initialize_shm_vec(n, DATA); 
                    sigsuspend(&mask);

                    sum = child_sum(vector, indexes[j*2], indexes[j*2+1]);
                    results[i] = sum;

                    shmdt(indexes);
                    shmdt(results);
                    shmdt(vector);
                    return EXIT_SUCCESS;
                default:
                    kids[j] = pid;
                    if((j+1) == subpr_n) {
                        struct timespec ts;
                        ts.tv_sec = 0;
                        ts.tv_nsec = 1e6;
                        nanosleep(&ts, &ts);
                    }
            }
        }

        vector = read_from_file(f, n);

        fill_indexes(n, subpr_n, indexes);

        for(int j = 0; j < subpr_n; j++) {
            results[j] = 0;
            kill(kids[i], SIGUSR1);
        }

        for (int j = 0; j < subpr_n; j++){
			wait(0);
		}

        for (int j = 0; j < subpr_n; j++){
			final_sum += results[j];
		}

        timer = clock() - timer;
        double time = (double) timer / CLOCKS_PER_SEC;
        times[i] = time;

		key_t key_idx = ftok(PROG_NAME, KEY_IDX);
		int size_idx = sizeof(int) * subpr_n * 2;
		int seg_idx = shmget(key_idx, size_idx, 0666 | IPC_CREAT);

		key_t key_res = ftok(PROG_NAME, KEY_RES);
		int size_res = sizeof(double) * subpr_n;
		int seg_res = shmget(key_res, size_res, 0666 | IPC_CREAT);

		key_t key_vec = ftok(PROG_NAME, KEY_DAT);
		int size_vec = sizeof(double) * n;
		int seg_vec = shmget(key_vec, size_vec, 0666 | IPC_CREAT);

		printf("Sum: %lf\n", final_sum);
		shmdt(indexes);
		shmdt(results);
		shmdt(vector);
		shmctl(seg_idx, IPC_RMID, NULL);
		shmctl(seg_res, IPC_RMID, NULL);
		shmctl(seg_vec, IPC_RMID, NULL);
        free(kids);
    }

    printf("~~~For vec with len: %d\n", n);
	for (int i = 0; i < SUBPR_ARR_LEN; i++)
		printf("~~~Calc time: %f, Subprocesses amount%d\n", times[i], subpr_arr[i]);

	free(times);
	return 0;

}