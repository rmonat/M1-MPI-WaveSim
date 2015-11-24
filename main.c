// TODO : iteration et steps plutôt des unsigned int
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <getopt.h>

#include "automaton.h"

void print_usage()
{
    printf("Usage: \n\n");
    
// TODO : et frprintf sur stderr en plus

}


typedef struct inst
{
    int iteration;
    double dt;
    char *input_path;
    char *lastdump;
    char *alldump;
    char *sensors;
    int x;
    int y;
}
inst;

void step0(inst instance, grid *g)
{
    int n = g->n;
    int m = g->m;
    int u;
    int v;
    double sqspeed = (g->v)*(g->v);
    grid tmp, current;
    new_grid(&tmp, g->n, g->m, g->v);
    new_grid(&current, g->n, g->m, g->v);

    copy_grid(g, &current);
    for(int s = 0; s < instance.iteration; s++)
    {
	for(size_t i = 0; i < g->n; i++)
	{
	    for(size_t j = 0; j < g->m; j++)
	    {
		u = current.data[j+i*n].u;
		v = current.data[j+i*n].v;
		tmp.data[j+i*n].u = u + v;
		tmp.data[j+i*n].v = v + sqspeed * (current.data[j+((i+1) % n)*n].u + current.data[j+((i-1) % n)*n].u + current.data[((j+1) % m) + i*n].u + current.data[((j-1) % m) + i*n].u - 4u);
	    }
	}

	copy_grid(&tmp, &current);
    }

    free(tmp.data);
    free(current.data);
	
}


int main(int argc, char** argv)
{
    int opt;
    int step = -1;
    int iteration = -1;
    double dt = 0;
    char *input_path = NULL;
    char *lastdump = NULL;
    char *alldump = NULL;
    char *sensor = NULL;
    int x = -1;
    int y = -1;
    
    static struct option long_options[] = {
        {"step",       required_argument, 0,  's' },
        {"i",          required_argument, 0,  'i' },
        {"iterations", required_argument, 0,  't' },
        {"dt",         required_argument, 0,  'd' },
	{"grid",       required_argument, 0,  'g' },
	{"lastdump",   required_argument, 0,  'l' },
	{"alldump",    required_argument, 0,  'a' },
	{"sensor",     required_argument, 0,  'S' },
        {0,           0,                 0,  0   }
    };


    int long_index = 0;
    while ((opt = getopt_long_only(argc, argv,"s:i:t:d:g:l:a:S:",
			      long_options, &long_index )) != -1) {
        switch (opt) {
	case 's' :
	    step = atoi(optarg);
	    break;
	case 'i' :
	    input_path = optarg;
	    break;
	case 't' :
	    iteration = atoi(optarg); 
	    break;
	case 'd' :
	    dt = (double) atof(optarg);
	    break;
	case 'g' :
	    x = atoi(optarg);
//	    printf("%c %c %d\n", argv[optind][0],argv[optind][1], strlen(argv[optind]));
	    if(optind < argc && argv[optind][0] != '-' && (strlen(argv[optind]) < 2 || argv[optind][1] != '-'))
	    {
		y = atoi(argv[optind]);
		optind++;
	    }
	    else
	    {
		fprintf(stderr, "Error: grid need to positive integers\n");
		print_usage();
	    }
	    break;
	case 'l':
	    lastdump = optarg;
	    break;
	case 'a':
	    alldump = optarg;
	    break;
	case 'S':
	    sensor = optarg;
	    break;
	default:
	    print_usage(); 
	    exit(EXIT_FAILURE);
        }
    }

#ifdef DEBUG
    fprintf(stderr, "Option Recap:\n"
	    "\tstep: %d\n"
	    "\tinput: %s\n"
	    "\titeration: %d\n"
	    "\tdt: %f\n"
	    "\tgrid: %d %d\n",
	    step, input_path, iteration, dt, x, y);
    if(lastdump != NULL)
    {
	fprintf(stderr, "\tlastdump: %s\n", lastdump);
//	char *b;
//	sprintf(b, lastdump, 2);
//	printf("\t\t%s\n", b);
    }
    if(alldump != NULL)
	fprintf(stderr, "\talldump: %s\n", alldump);
    if(sensor != NULL)
	fprintf(stderr, "\tsensor: %s\n", sensor);

#endif

    int err = 0;
    // Check we have at least the required arguments.
    if(!(step >= 0 && step <= 4))
    {
	err = 1;
	fprintf(stderr, "Error: step must be between 0 and 4\n");
    }
    if(iteration <= 0)
    {
	err = 1;
	fprintf(stderr, "Error: the number of iterations should be a positive integer\n");
    }
    if(dt == 0)
    {
	err = 1;
	fprintf(stderr, "Error: dt should not be 0\n");
    }
    if(x <= 0)
    {
	err = 1;
	fprintf(stderr, "Error: the x-coordinate of the grid should be a positive integer\n");
    }
    if(y <= 0)
    {
	err = 1;
	fprintf(stderr, "Error: the y-coordinate of the grid should be a positive integer\n");
    }

    if(err == 1)
    {
	print_usage(); 
	exit(EXIT_FAILURE);
    }

    // TESTS
    grid g;
    parse_file(input_path, &g);
    g.data[0].u = 1;
    g.data[1].u = 2;
    dump_grid("bla", &g);
    
    return 0;
}
