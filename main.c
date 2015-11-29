// TODO : iteration et steps plutôt des unsigned int
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <mpi.h>

#include "automaton.h"
#include "args.h"
#include "step0.h"


int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // TODO : gérer, n'utiliser qu'un rang
    inst instance;
    parse_options(&instance, argc, argv);

    if(instance.step == 0 && rank == 0)
    {
	grid g;
	parse_file(instance.input_path, &g);
	step0(instance, &g);
    }
    else
    {
	if(size != instance.p * instance.q && rank == 0)
	{
	    fprintf(stderr, "Error: there should be %d procs", (instance.p * instance.q));
	    print_usage();
	    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	}

	// Creation of the 2D torus we will then use
	MPI_Comm comm;
	int dim[2] = {instance.p, instance.q};
	int period[2] = {1, 1};
	int reorder = 1;
	int coord[2];
	MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &comm);
	MPI_Cart_coords(comm, rank, 2, coord);


	grid global_grid;
	char type = 0;

	MPI_File input_file;

	// We start by reading the header of the file
	MPI_File_open(comm,  instance.input_path, MPI_MODE_RDONLY, MPI_INFO_NULL, &input_file);
	MPI_File_read_all(input_file, &type, 1, MPI_CHAR, MPI_STATUS_IGNORE);
	// maybe we need to swap the next 2 lines
	MPI_File_read_all(input_file, &(global_grid.m), 1, MPI_UINT64_T, MPI_STATUS_IGNORE);
	MPI_File_read_all(input_file, &(global_grid.n), 1, MPI_UINT64_T, MPI_STATUS_IGNORE);
	MPI_File_read_all(input_file, &(global_grid.v), 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

#ifdef DEBUG
	if(rank == 0)
	    printf("n, m, v = %d %d %f\n", global_grid.n, global_grid.m, global_grid.v);
#endif

	size_t nb_cells = global_grid.n * global_grid.m;

	// Now we create the data structures.
	int blocks[2] = {1, 1};
	MPI_Datatype types[2] = {MPI_BYTE, MPI_DOUBLE};
	MPI_Aint p_size = 9;
	MPI_Aint p_disp[2] = {0, 1};
	MPI_Aint a_size = sizeof(cell);
	MPI_Aint a_disp[2] = {offsetof(cell, type), offsetof(cell, u)};

	MPI_Datatype p_tmp, a_tmp, p_cell, a_cell;

	// Aligned struct, memory representation
	MPI_Type_create_struct(2, blocks, a_disp, types, &a_tmp);
	MPI_Type_create_resized(a_tmp, 0, a_size, &a_cell);
	MPI_Type_commit(&a_cell);
	    
	// Packed struct, file-based representation
	MPI_Type_create_struct(2, blocks, p_disp, types, &p_tmp);
	MPI_Type_create_resized(p_tmp, 0, p_size, &p_cell);
	MPI_Type_commit(&p_cell);
	

	// Now, we create our matrix
	MPI_Datatype matrix;
	int sizes[2] = {global_grid.n, global_grid.m};
	int subsizes[2] = {global_grid.n / instance.p, global_grid.m / instance.q};
	int starts[2] = {0, 0};
	MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, p_cell, &matrix);
	MPI_Type_commit(&matrix);

	// Set file view for each element
	MPI_Offset grid_start;
	MPI_File_get_position(input_file, &grid_start);


	size_t local_size = nb_cells / (instance.p * instance.q);
	
	MPI_File_set_view(input_file, grid_start + global_grid.m*global_grid.n/(instance.p)*p_size*coord[0] + global_grid.m/(instance.q)*p_size*coord[1], p_cell, matrix, "native", MPI_INFO_NULL);

	// allocate the cell array we will use
	cell *cells = malloc(local_size*sizeof(cell));

	MPI_File_read_all(input_file, cells, local_size, a_cell, MPI_STATUS_IGNORE);


#ifdef DEBUG
	for(size_t i = 0; i < local_size; i++)
	    fprintf(stderr, "%d - %d %f\n", rank, cells[i].type, cells[i].u);
#endif


	// Some cleaning
	free(cells);
	MPI_File_close(&input_file);
	MPI_Type_free(&p_cell);
	MPI_Type_free(&a_cell);
	MPI_Type_free(&matrix);
	MPI_Finalize();
    }
    return 0;
}

