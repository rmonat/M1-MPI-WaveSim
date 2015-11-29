typedef struct inst
{
    int iteration;
    int step;
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
void parse_options(inst *instance, int argc, char **argv);