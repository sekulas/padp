#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define N 10000000
#define S (int)sqrt(N)
#define M N/10
#define THREADS_N 8
#define CHUNK_SIZE 1000

int main(int argc, char**argv){
    long int a[S + 1]; /*helping array*/
    long int primes[M]; /*prime numbers 2..N*/
    
    long i, k, number, rest;
    long int aliquot_n; /*n of aliquots*/
    long int primes_n = 0; /*prime numbers amount in primes arr*/
    double time;
    FILE *fp;

    omp_set_num_threads(THREADS_N);
    int chunk_size = S / THREADS_N;

    /*looking for aliquots in 2..S*/
    #pragma omp parallel for schedule(dynamic, chunk_size) private(i)
    for(i = 2; i <= S; i++) {
        a[i] = 1;
    }

    for(i=2; i<=S; i++) {
        if(a[i] == 1){
            primes[primes_n++] = i; /*aliquot saving*/
            /*removing numbers being multiplication of i*/
            
            #pragma omp parallel for schedule(dynamic, chunk_size) private(k)
            for(k = i+i; k<=S; k+=i) {
                a[k] = 0;
            }
        }
    }

	chunk_size = (int)(N-S+1)/(CHUNK_SIZE);
    aliquot_n = primes_n; /*remembering the amount of aliquotes*/

    /*finding primes*/
    #pragma omp parallel for schedule(dynamic, chunk_size) private(number, rest, k)
    for(number = S+1; number <= N; number++){
        for(k = 0; k < aliquot_n; k++){
            rest = (number % primes[k]);
            if(rest == 0) {
                break; /*composite number*/
            }
        }
        if(rest != 0) {
            #pragma omp critical
            {
                primes[primes_n++] = number; /*saving prime number*/
            }
        }
    }

    if((fp = fopen("primes.txt", "w")) == NULL) {
        printf("Cannot open the file for writing.\n");
        exit(1);
    }

    for(i=0; i< primes_n; i++) {
        fprintf(fp,"%ld ", primes[i]);
    }

    fclose(fp);
    return 0;
}