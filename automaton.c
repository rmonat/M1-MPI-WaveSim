#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "automaton.h"

// TODO regarder les return values des fwrite, fread ?

void new_grid(grid* g, size_t n, size_t m, double v)
{
    g->n = n;
    g->m = m;
    g->v = v;
    g->data = (cell *) calloc (n*m, sizeof(cell));
}

void copy_grid(grid* src, grid* dst)
{
    assert(src->n == dst->n && src->m == dst->m && src->v == dst->v);
    memcpy(dst->data, src->data, sizeof(double)*(src->n)*(src->m));
}

void dump_grid(char* filename, grid *g)
{
    FILE *fp;
    fp = fopen(filename, "wb");
    if(fp == NULL)
    {
	perror("Error during grid dumping:");
    }

    int n = g->n;
    for(size_t i = 0; i < g->n; i++)
    {
	for(size_t j = 0; j < g->m; j++)
	{
	    fwrite(&(g->data[i*n+j].u), sizeof(double), 1, fp);
	    // TODO : danger, gérer les walls
	}
    }

    fclose(fp);

}

void parse_file(char* filename, grid *g)
{
    // TODO : ET LES AUTRES STEP ?!

    FILE *fp;
    fp = fopen(filename, "rb");
    if(fp == NULL)
    {
	perror("Error during the input parsing: ");
    }
    
    char test;
    size_t n;
    size_t m;
    double v;
    fread(&test, sizeof(char), 1, fp);
    fread(&n, sizeof(size_t), 1, fp);
    fread(&m, sizeof(size_t), 1, fp);
    fread(&v, sizeof(double), 1, fp);

    char cell_type;
    double cell_value;

    g->n = n;
    g->m = m;
    g->v = v;
    g->data = (cell *) calloc(n*m, sizeof(cell));
    
    for(size_t i = 0; i < n; i++)
	for(size_t j = 0; j < m; j++)
	{
	    fread(&cell_type, sizeof(char), 1, fp);
	    fread(&cell_value, sizeof(double), 1, fp);
	    g->data[i*n+j].type = 1 - cell_type; // my definition, cells are 1, walls are 0
	    g->data[i*n+j].u = cell_value;
	    g->data[i*n+j].v = 0;
	}

#ifdef DEBUG
    printf("Initial grid: %d %zd %zd %f\n", test, n, m, v);
    for(size_t i = 0; i < n; i++)
    {
	for(size_t j = 0; j < m; j++)
	    printf("%f ", g->data[i*n+j].u);
	printf("\n");
    }
#endif
    
    fclose(fp);
}
