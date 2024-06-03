EXECS=v_sum_scat
MPICC?=mpicc
MPIEXEC?=mpiexec

gen: gen_vec.c
	${MPICC} -o gen_vec gen_vec.c

gen_ex:
	./gen_vec

scat: v_sum_scat.c
	${MPICC} -o v_sum_scat v_sum_scat.c

scat_ex:
	${MPIEXEC} -n 5 ./v_sum_scat

xd: v_sum_scat2.c
	${MPICC} -o v_sum_scat2 v_sum_scat2.c

xd_ex:
	${MPIEXEC} -n 5 ./v_sum_scat2