#ifndef ARGS
#define ARGS

typedef struct inst
{
    int iteration;
    int step;
    int frequency;
    double dt;
    char *input_path;
    char *lastdump;
    char *alldump;
    char *sensors;
    int p;
    int q;
}
inst;


void print_usage();
void init_inst(inst *instance);
int parse_options(inst *instance, int argc, char **argv, int rank);

#endif
