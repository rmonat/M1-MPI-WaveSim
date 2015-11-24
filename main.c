#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <getopt.h>

void print_usage()
{
    printf("Usage: \n\n");
    
// TODO : et frprintf sur stderr en plus
}

int main(int argc, char** argv)
{
    int opt;
    int step = 0;
    int iteration = 0;
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
    while ((opt = getopt_long(argc, argv,"s:i:t:d:g:l:a:S:",
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
	fprintf(stderr, "\tlastdump: %s\n", lastdump);
    if(alldump != NULL)
	fprintf(stderr, "\talldump: %s\n", alldump);
    if(sensor != NULL)
	fprintf(stderr, "\tsensor: %s\n", sensor);

#endif

    // Check we have at least the required arguments.
    return 0;
}
