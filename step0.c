#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "automaton.h"
#include "args.h"

void step0(inst instance, grid *g)
{
    size_t n = g->n;
    size_t m = g->m;
    char *b = malloc(256);
    double sqspeed = (g->v)*(g->v);
    grid tmp, current;
    new_grid(&tmp, g->n, g->m, g->v);
    new_grid(&current, g->n, g->m, g->v);

#ifdef DEBUG
    printf("sqspeed: %f\nn, m = %zu %zu\n", sqspeed, n, m);
#endif
    
    copy_grid(g, &current);

    dump_grid("current", &current);
    printf("\t\t\t%f\n", current.data[34+39*256].u);
    
    for(int s = 0; s < instance.iteration; s++)
    {
	for(size_t i = 0; i < n; i++)
	{
	    for(size_t j = 0; j < m; j++)
	    {
		// TODO : gérer les bords
		tmp.data[j+i*m].u = current.data[j+i*m].u + (current.data[j+i*m].v * instance.dt);
		tmp.data[j+i*m].v = current.data[j+i*m].v + sqspeed * (current.data[j+((i+1) % n)*m].u + current.data[j+((n+i-1) % n)*m].u + current.data[((j+1) % m) + i*m].u + current.data[((m+j-1) % m) + i*m].u - (4 * current.data[j+i*m].u)) * instance.dt;
		
	    }
	}
//	printf("\t\t\t%f\n", current.data[158+42*256].u);

	if(instance.alldump != NULL)
	{
	    sprintf(b, instance.alldump, s);
	    dump_grid(b, &tmp);
	}
//	printf("%f %f -> %f %f\n", current.data[38+38*n].u, current.data[38+38*n].v, tmp.data[38+38*n].u, tmp.data[38+38*n].v);
	copy_grid(&tmp, &current);
    }

    if(instance.lastdump != NULL)
    {
	dump_grid(instance.lastdump, &current);
    }
    
    free(tmp.data);
    free(current.data);
    free(b);
}
