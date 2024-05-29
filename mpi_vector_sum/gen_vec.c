#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char** argv) {
    int count = argc > 1 ? atoi(argv[1]) : 10;
    double lower = argc > 2 ? atof(argv[2]) : 0;
    double upper = argc > 3 ? atof(argv[3]) : 1;
    FILE *out = argc > 4 ? fopen(argv[4], "w") : fopen("vector.dat", "w");

    int i;
    double sum = 0.0f;

    if(out == NULL) {
        fprintf(stderr, "I cannot write to a file.\n");
    }

    srand(time(0));
    fprintf(out, "%d\n", count);

    for(i = 0; i < count; i++) {
        double num = ((double)rand() / RAND_MAX) * (upper - lower) + lower;
        sum += num;
        fprintf(out, "%f\n", num);
    }

    fclose(out);

    printf("Generated vector elems sum: %f\n", sum);
}