#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <assert.h>
#include <stddef.h>

int main(int argc, char** argv)
{
    typedef struct
    {
	char type;
	double val;
    } input;

    int rank, size;
    
    MPI_Init(&argc, &argv);
    MPI_Comm comm;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int dim[2] = {2, 2};
    int period[2] = {0, 0};
    int reorder = 1;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &comm);
    int coord[2];
    MPI_Cart_coords(comm, rank, 2, coord);

    assert(size == 4);

    MPI_File in;
    MPI_Offset filesize;
    MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &in);
    MPI_File_get_size(in, &filesize);

    int mysize = filesize/size;
    if(rank == 0)
	printf("filesize, mysize %d %d\n", filesize, mysize);
    int globalstart = rank * mysize;

    
    
    int blocks[2] = {1,1};
    MPI_Datatype types[2] = {MPI_BYTE, MPI_DOUBLE};
    MPI_Aint displacements[2];
    MPI_Datatype cell_type, cell_packed, cp;
    displacements[0] = offsetof(input, type);
    displacements[1] = offsetof(input, val);
    MPI_Aint packed_disp[2];
    packed_disp[0] = 0;
    packed_disp[1] = 1;
    MPI_Type_create_struct(2, blocks, displacements, types, &cell_type);
    MPI_Type_commit(&cell_type);

    MPI_Type_create_struct(2, blocks, packed_disp, types, &cp);
    MPI_Type_create_resized(cp, 0, 9, &cell_packed);
    MPI_Type_commit(&cell_packed);
    int cell_type_size;
    MPI_Type_size(cell_type, &cell_type_size);
    input chunk[4];// = malloc(4*sizeof(input));//cell_type_size);


    MPI_Datatype matrix;
    int sizes[2] = {4, 4};
    int subsizes[2] = {2, 2};
    int starts[2] = {0, 0};

    MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, cell_type, &matrix);
    MPI_Type_commit(&matrix);
    int matrix_type_size;
    MPI_Type_size(matrix, &matrix_type_size);
    if(rank == 0)
    	printf("matsize %d, bla %d\n", matrix_type_size, mysize/cell_type_size);

    MPI_Datatype test;
    MPI_Type_contiguous(4, cell_packed, &test);
    MPI_Type_commit(&test);
    
//    MPI_File_set_view(in, 4*rank, cell_type, cell_type, "native", MPI_INFO_NULL);
//    MPI_File_set_view(in, 0, cell_packed, matrix, "native", MPI_INFO_NULL);
    MPI_File_set_view(in, 0, cell_packed, test, "native", MPI_INFO_NULL);
    MPI_File_read_at_all(in, 0, chunk, 4, cell_type, MPI_STATUS_IGNORE);

    
    printf("%d - Got %f %d %f %d %f %d %f %d\n", rank, chunk[0].val, chunk[0].type);//, chunk[1].val, chunk[1].type, chunk[2].val, chunk[2].type, chunk[3].val, chunk[3].type) ;
    
    MPI_File_close(&in);
    MPI_Finalize();
}
