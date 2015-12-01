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
    if(parse_options(&instance, argc, argv, rank) == 1)
    {
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	exit(EXIT_FAILURE);
    }
//	MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

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
	// we needed to swap the next 2 lines
	MPI_File_read_all(input_file, &(global_grid.n), 1, MPI_UINT64_T, MPI_STATUS_IGNORE);
	MPI_File_read_all(input_file, &(global_grid.m), 1, MPI_UINT64_T, MPI_STATUS_IGNORE);
	MPI_File_read_all(input_file, &(global_grid.v), 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

#ifdef DEBUG
	if(rank == 0)
	    printf("n, m, v = %d %d %f\n", global_grid.n, global_grid.m, global_grid.v);
#endif

	size_t nb_cells = global_grid.n * global_grid.m;

	// Now we create the data structures.
	int blocks[2] = {1, 1};
	MPI_Datatype types[2] = {MPI_BYTE, MPI_DOUBLE};
	MPI_Aint a_size = sizeof(cell);
	MPI_Aint a_disp[2] = {offsetof(cell, type), offsetof(cell, u)};

	MPI_Aint p_size = 9;
	MPI_Aint p_disp[2] = {0, 1};

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

	// We extend this matrix
	MPI_Datatype ematrix;
	int e_subsizes[2] = {2 + subsizes[0], 2 + subsizes[1]};
	int e_start[2] = {1, 1};
	MPI_Type_create_subarray(2, e_subsizes, subsizes, e_start, MPI_ORDER_C, a_cell, &ematrix);
	MPI_Type_commit(&ematrix);
	

	// The next 3 types are for the export of the grid
	MPI_Datatype d_type;
	MPI_Type_create_resized(MPI_DOUBLE, 0, sizeof(cell), &d_type);
	MPI_Type_commit(&d_type);
	

	MPI_Datatype d_matrix;
	MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &d_matrix);
	MPI_Type_commit(&d_matrix);

	MPI_Datatype d_rmatrix; // to go from the extended matrix with ghost zones to the other one
	MPI_Type_create_subarray(2, e_subsizes, subsizes, e_start, MPI_ORDER_C, d_type, &d_rmatrix);
	MPI_Type_commit(&d_rmatrix);



	// Set file view for each element
	MPI_Offset grid_start;
	MPI_File_get_position(input_file, &grid_start);


	size_t local_size = nb_cells / (instance.p * instance.q);
	size_t local_nrows = global_grid.n/instance.p;
	size_t local_ncols = global_grid.m/instance.q;
	
	MPI_File_set_view(input_file, grid_start + global_grid.m*global_grid.n/(instance.p)*p_size*coord[0] + global_grid.m/(instance.q)*p_size*coord[1], p_cell, matrix, "native", MPI_INFO_NULL);

	// allocate the cell array we will use
	cell **cells;
	cells = malloc(2*sizeof(cell *));

	
	cells[1] = malloc((2+global_grid.n/instance.p)*(2+global_grid.m/instance.q)*sizeof(cell));
	cells[0] = malloc((2+global_grid.n/instance.p)*(2+global_grid.m/instance.q)*sizeof(cell));

	MPI_File_read_all(input_file, cells[0], 1, ematrix, MPI_STATUS_IGNORE);



#ifdef DEBUG
	for(size_t i = 1; i < 1+global_grid.n/instance.p; i++)
	    for(size_t j = 1; j < 1+global_grid.m/instance.q; j++)
		fprintf(stderr, "%d - %d %f\n", rank, cells[0][i*(2+global_grid.m/instance.q)+j].type, cells[0][i*(2+global_grid.m/instance.q)+j].u);
#endif

	MPI_Datatype l_row; // local row
	MPI_Type_contiguous(global_grid.n/instance.p, a_cell, &l_row);
	MPI_Type_commit(&l_row);

	MPI_Datatype l_col; // local column. A bit trickier, we need a type_vector.
	MPI_Type_vector(local_nrows, 1, local_ncols+1, a_cell, &l_col);
	MPI_Type_commit(&l_col);

	
	int top, bot, left, right, t[2], b[2], r[2], l[2];
	double sqspeed = global_grid.v * global_grid.v;

	int curr = 0, next = 0;
	char *alldump = malloc(256);

	if(rank == 0)
	{
	    /* cells[0][2+local_ncols+1].u = 1; */
	    /* cells[0][2+local_ncols+2].u = 1; */
	    /* cells[0][2+local_ncols+3].u = 1;  */
	    /* cells[0][2+local_ncols+4].u = 1;  */
	    /* cells[0][2*(2+local_ncols)+1].u = 1; */
	    /* cells[0][2*(2+local_ncols)+2].u = 1; */
	    /* cells[0][2*(2+local_ncols)+3].u = 1; */
//	    cells[0][2*(2+local_ncols)+4].u = 1;

	}


	for(int s = 0; s < instance.iteration; s++)
	{
	    //  TODO : SYNCHRONISER MIEUX ? Ou peut être même pas besoin de synchroniser en fait
	    MPI_Barrier(MPI_COMM_WORLD);

	    // on va communiquer dans cell[curr]...
	    // et mettre à jour dans cell[next]...
	    curr = s % 2;
	    next = (s+1) % 2;
	    
	    // We copy the edges of the grid.
	    // We first need the ranks of the neighbours
	    t[0] = (coord[0]+1) % instance.p;
	    t[1] = coord[1];
	
	    b[0] = (instance.p+coord[0]-1) % instance.p;
	    b[1] = coord[1];

	    r[0] = coord[0];
	    r[1] = (coord[1]+1) % instance.q;
	    
	    l[0] = coord[0];
	    l[1] = (instance.q+coord[1]-1) % instance.q;

	    MPI_Cart_rank(comm, t, &top);
	    MPI_Cart_rank(comm, b, &bot);
	    MPI_Cart_rank(comm, r, &right);
	    MPI_Cart_rank(comm, l, &left);

	    // Then we need to update the edges of our local grid
	    // Update top and bottom rows
	    MPI_Sendrecv(&(cells[curr][local_ncols+3]),                     1, l_row, top, 0,
			 &(cells[curr][(local_ncols+2)*(local_nrows+1)+1]), 1, l_row, bot, 0,
			 comm, MPI_STATUS_IGNORE);
	
	    MPI_Sendrecv(&(cells[curr][(local_ncols+2)*(local_nrows)+1]),   1, l_row, bot, 0,
			 &(cells[curr][1]),                                 1, l_row, top, 0,
			 comm, MPI_STATUS_IGNORE);
	
	    // Update left and right
	    MPI_Sendrecv(&(cells[curr][local_ncols+3]),                     1, l_col, left,  0,
			 &(cells[curr][2*(local_ncols+2)-1]),               1, l_col, right, 0,
			 comm, MPI_STATUS_IGNORE);

	    MPI_Sendrecv(&(cells[curr][2*(local_ncols+2)-2]),               1, l_col, right, 0,
			 &(cells[curr][local_ncols+2]),                     1, l_col, left,  0,
			 comm, MPI_STATUS_IGNORE);



	    // We compute the update of the grid
	    for(size_t i = 1; i < 1+local_nrows; i++)
	    {
		for(size_t j = 1; j < 1+local_ncols; j++)
		{
		    cells[next][j+i*(2+local_ncols)].u = cells[curr][j+i*(2+local_ncols)].u + (cells[curr][j+i*(2+local_ncols)].v * instance.dt);
		    cells[next][j+i*(2+local_ncols)].v = cells[curr][j+i*(2+local_ncols)].v + sqspeed * (cells[curr][j+(i+1)*(2+local_ncols)].u + cells[curr][j+(i-1)*(2+local_ncols)].u + cells[curr][(j+1) + i*(2+local_ncols)].u + cells[curr][(j-1) + i*(2+local_ncols)].u - (4 * cells[curr][j+i*(2+local_ncols)].u)) * instance.dt;
		    
		}
	    }

	    if(instance.alldump != NULL)
	    {
		MPI_File dump_file;

		sprintf(alldump, instance.alldump, s);
		MPI_File_open(comm, alldump, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &dump_file);
		
		MPI_File_set_view(dump_file, global_grid.m*global_grid.n/(instance.p)*sizeof(double)*coord[0] + global_grid.m/(instance.q)*sizeof(double)*coord[1], MPI_DOUBLE, d_matrix, "native", MPI_INFO_NULL);
		
		MPI_File_write_all(dump_file, &(cells[curr][0].u), 1, d_rmatrix, MPI_STATUS_IGNORE);
		MPI_File_close(&dump_file);


	    }
	}

	
	if(instance.lastdump != NULL)
	{
	    // bon, comment on fait ça ? peut être qu'en faisant un resize ça marche ?
	    MPI_File last_file;
	    MPI_File_open(comm, instance.lastdump, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &last_file); 
	    MPI_File_set_view(last_file, global_grid.m*global_grid.n/(instance.p)*sizeof(double)*coord[0] + global_grid.m/(instance.q)*sizeof(double)*coord[1], MPI_DOUBLE, d_matrix, "native", MPI_INFO_NULL); // déjà, il y a un grid_strat en trop, d_type ou MPI_DOUBLE ?

	    MPI_File_write_all(last_file, &(cells[next][0].u), 1, d_rmatrix, MPI_STATUS_IGNORE);
	    MPI_File_close(&last_file);


	}
	

	// Some cleaning
	free(cells);
	free(alldump);
	MPI_File_close(&input_file); // TODO : mettre plus tôt
	MPI_Type_free(&p_cell);
	MPI_Type_free(&a_cell);
	MPI_Type_free(&matrix);
	MPI_Finalize();
    }
    return 0;
}

// on fait le truc de la non copie et alternance modulo 2 ?
//Sendrecv des edges
//void update_grid(cell *cells, quelques éléments de l'instance);

