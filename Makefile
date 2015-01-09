mozco_no_espacial: ciclo_mozco_no_espacial.c libPP_5.0.c EntSalArb_MP.c GNA.c MC_sweep_mozquito.c
	gcc -o3 -fopenmp ciclo_mozco_no_espacial.c libPP_5.0.c EntSalArb_MP.c GNA.c MC_sweep_mozquito.c -lfftw3 -lm -o mozco.out
