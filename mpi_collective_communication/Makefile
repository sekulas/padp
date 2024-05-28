EXECS=gcd pipeline
MPICC?=mpicc
MPIEXEC?=mpiexec

gcd: gcd.c
	${MPICC} -o gcd gcd.c -lm

gcd_ex:
	${MPIEXEC} -n 5 ./gcd 3 6 9 12 15

pipeline: pipeline.c
	${MPICC} -o pipeline pipeline.c

pipeline_ex:
	${MPIEXEC} -n 5 ./pipeline

clean:
	rm -f ${EXECS}
