#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>

typedef struct {
    double u;
    char type;
    double v;
} cell;

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 3) {
        fprintf(stderr, "Usage: %s p q", argv[0]);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    int p = atoi(argv[1]);
    int q = atoi(argv[2]);
    assert(size == p*q);
    
    MPI_Comm comm;
    int dim[2] = {p, q};
    int period[2] = {1, 1};
    int reorder = 1;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &comm);
    int coord[2];
    MPI_Cart_coords(comm, rank, 2, coord);

    size_t height = 4;
    size_t width = 4;
    cell *cells = malloc(height/p*width/q*sizeof(cell));

    for(int i = 0; i < height*width/(p*q); i++)
    {
	cells[i].type = 1;
	cells[i].u = rank + 0.5;
	cells[i].v = width*height - rank;
    }

    MPI_Datatype d_type;
    MPI_Type_create_resized(MPI_DOUBLE, 0, sizeof(cell), &d_type);
    MPI_Type_commit(&d_type);

    MPI_Datatype mem_matrix;
    MPI_Type_contiguous(width*height/(p*q), d_type, &mem_matrix);
    MPI_Type_commit(&mem_matrix);

    int sizes[2] = {height, width};
    int subsizes[2] = {height/p, width/q};
    int starts[2] = {0, 0};
   
    MPI_Datatype d_matrix;
    MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &d_matrix);
    MPI_Type_commit(&d_matrix);



    MPI_File last_file;
    char *lastdump = "bla.dump";
    MPI_File_open(comm, lastdump, MPI_MODE_WRONLY | MPI_MODE_CREATE , MPI_INFO_NULL, &last_file); 

    MPI_File_set_view(last_file, width*height/p*sizeof(double)*coord[0]+width/q*sizeof(double)*coord[1], MPI_DOUBLE, d_matrix, "native", MPI_INFO_NULL);

    MPI_File_write_all(last_file, &(cells[0]), 1, mem_matrix, MPI_STATUS_IGNORE);

   
    free(cells);
    MPI_Type_free(&mem_matrix);
    MPI_Type_free(&d_matrix);
    MPI_Type_free(&d_type);
    MPI_File_close(&last_file);
    MPI_Finalize();
    return 0;
}
