// TODO : iteration et steps plutôt des unsigned int
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <mpi.h>

#include "automaton.h"
#include "args.h"
#include "step0.h"
#include "step123.h"

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    inst instance;
    if(parse_options(&instance, argc, argv, rank) == 1)
    {
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	exit(EXIT_FAILURE);
    }

    if(instance.step == 0)
    {
	if(rank == 0)
	{
	    grid g;
	    parse_file(instance.input_path, &g);
	    step0(instance, &g);
	}
    }
    else
    {
	if(size != instance.p * instance.q && rank == 0)
	{
	    fprintf(stderr, "Error: there should be %d procs", (instance.p * instance.q));
	    print_usage();
	    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	}

	if(instance.step <= 3)
	    step123(instance, rank, size);
	else if(instance.step == 4)
	    step4(instance, rank, size);
	
	MPI_Finalize();
    }
    return 0;
}

// on fait le truc de la non copie et alternance modulo 2 ?
//Sendrecv des edges
//void update_grid(cell *cells, quelques éléments de l'instance);

