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
    //grid g;
    //parse_file(input_path, &g);
    //g.data[0].u = 1;
    //g.data[1].u = 2;
    //dump_grid("bla", &g);


    if(step == 0)
    {
	inst instance;
	instance.iteration = iteration;
	instance.dt = dt;
	instance.input_path = input_path;
	instance.lastdump = lastdump;
	instance.alldump = alldump;
	instance.sensors = sensor;
	instance.x = x;
	instance.y = y;

	grid g;
	parse_file(input_path, &g);
	step0(instance, &g);
    }
    return 0;
}

