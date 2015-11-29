#include <getopt.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "args.h"

void print_usage()
{
    printf("Usage: \n\n");
    
// TODO : et frprintf sur stderr en plus

}

void init_inst(inst *instance)
{
    instance->iteration = -1;
    instance->step = -1;
    instance->dt = 0;
    instance->input_path = NULL;
    instance->lastdump = NULL;
    instance->alldump = NULL;
    instance->sensors = NULL;
    instance->p = -1;
    instance->q = -1;
}

void parse_options(inst *instance, int argc, char **argv)
{
    int opt;
    init_inst(instance);
    
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
	    instance->step = atoi(optarg);
	    break;
	case 'i' :
	    instance->input_path = optarg;
	    break;
	case 't' :
	    instance->iteration = atoi(optarg); 
	    break;
	case 'd' :
	    instance->dt = (double) atof(optarg);
	    break;
	case 'g' :
	    instance->p = atoi(optarg);
//	    printf("%c %c %d\n", argv[optind][0],argv[optind][1], strlen(argv[optind]));
	    if(optind < argc && argv[optind][0] != '-' && (strlen(argv[optind]) < 2 || argv[optind][1] != '-'))
	    {
		instance->q = atoi(argv[optind]);
		optind++;
	    }
	    else
	    {
		fprintf(stderr, "Error: grid need to positive integers\n");
		print_usage();
		exit(EXIT_FAILURE);
	    }
	    break;
	case 'l':
	    instance->lastdump = optarg;
	    break;
	case 'a':
	    instance->alldump = optarg;
	    break;
	case 'S':
	    instance->sensors = optarg;
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
	    instance->step, instance->input_path, instance->iteration, instance->dt, instance->p, instance->q);
    if(instance->lastdump != NULL)
    {
	fprintf(stderr, "\tlastdump: %s\n", instance->lastdump);
//	char *b;
//	sprintf(b, lastdump, 2);
//	printf("\t\t%s\n", b);
    }
    if(instance->alldump != NULL)
	fprintf(stderr, "\talldump: %s\n", instance->alldump);
    if(instance->sensors != NULL)
	fprintf(stderr, "\tsensor: %s\n", instance->sensors);

#endif

    int err = 0;
    // Check we have at least the required arguments.
    if(!(instance->step >= 0 && instance->step <= 4))
    {
	err = 1;
	fprintf(stderr, "Error: step must be between 0 and 4\n");
    }
    if(instance->iteration <= 0)
    {
	err = 1;
	fprintf(stderr, "Error: the number of iterations should be a positive integer\n");
    }
    if(instance->dt == 0)
    {
	err = 1;
	fprintf(stderr, "Error: dt should not be 0\n");
    }
    if(instance->p <= 0)
    {
	err = 1;
	fprintf(stderr, "Error: the p-coordinate of the grid should be a positive integer\n");
    }
    if(instance->q <= 0)
    {
	err = 1;
	fprintf(stderr, "Error: the q-coordinate of the grid should be a positive integer\n");
    }

    if(err == 1)
    {
	print_usage(); 
	exit(EXIT_FAILURE);
    }
}
