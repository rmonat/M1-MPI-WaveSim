#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "args.h"
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
    memcpy(dst->data, src->data, sizeof(cell)*(src->n)*(src->m));
}

void dump_grid(char* filename, grid *g)
{
    FILE *fp;
    fp = fopen(filename, "wb");
    if(fp == NULL)
    {
	perror("Error during grid dumping:");
    }

    size_t n = g->n;
    size_t m = g->m;
    
    for(size_t i = 0; i < n; i++)
    {
	for(size_t j = 0; j < m; j++)
	{
	    fwrite(&(g->data[i*m+j].u), sizeof(double), 1, fp);
	}
    }

    fclose(fp);

}

void parse_file(char* filename, grid *g)
{
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
	    if(cell_type == 1)
		printf("%zu %zu\n", i, j);
	    g->data[i*m+j].type = cell_type;
	    g->data[i*m+j].u = cell_value;
	    g->data[i*m+j].v = 0;
	}

#ifdef DEBUG
    fprintf(stderr, "Initial grid: %d %zd %zd %f\n", test, n, m, v);
    for(size_t i = 0; i < n; i++)
    {
    	for(size_t j = 0; j < m; j++)
	{
    	    fprintf(stderr, "%f ", g->data[i*m+j].u);
	}
    	fprintf(stderr, "\n");
    }
#endif
    
    fclose(fp);
}
