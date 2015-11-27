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
	double val;
	char type;
    } input;

    int rank, size;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    assert(size == 4);

    MPI_File in;
    MPI_Offset filesize;
    MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &in);
    MPI_File_get_size(in, &filesize);

    int mysize = filesize/size;
    int globalstart = rank * mysize;

    input *chunk = malloc(sizeof(input)*2);
//    void *chunk = malloc(18);
    
    
    int blocks[2] = {1,1};
    MPI_Datatype types[2] = {MPI_BYTE, MPI_DOUBLE};
    MPI_Aint displacements[2];
    MPI_Datatype cell_type;
    displacements[0] = offsetof(input, type);
    displacements[1] = offsetof(input, val);
    MPI_Type_create_struct(2, blocks, displacements, types, &cell_type);
    MPI_Type_commit(&cell_type);
    int cell_type_size;
    MPI_Type_size(cell_type, &cell_type_size);

    MPI_File_read_at_all(in, globalstart, chunk, mysize / cell_type_size, cell_type, MPI_STATUS_IGNORE);
    printf("%d - Got %f %d %f %d\n", rank, chunk[0].val, chunk[0].type, chunk[1].val, chunk[1].type);
    
    MPI_File_close(&in);
    MPI_Finalize();
}
