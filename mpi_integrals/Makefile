EXECS=integrate_send integrate_isend integrate_scatter
MPICC?=mpicc
MPIEXEC?=mpiexec

send: integrate_send.c
	${MPICC} -o integrate_send integrate_send.c

send_ex:
	${MPIEXEC} -n 5 ./integrate_send 0 1 100

isend: integrate_isend.c
	${MPICC} -o integrate_isend integrate_isend.c

isend_ex:
	${MPIEXEC} -n 5 ./integrate_isend 0 1 100

scat: integrate_scatter.c
	${MPICC} -o integrate_scatter integrate_scatter.c

scat_ex:
	${MPIEXEC} -n 5 ./integrate_scatter 0 1 100


clean:
	rm -f ${EXECS}
