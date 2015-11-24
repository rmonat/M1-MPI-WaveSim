typedef struct cell
{
    char type;
    double u;
    double v; // v = u'
} cell;

typedef struct grid
{
    size_t n;
    size_t m;
    double v;
    cell *data;
}
grid;


void new_grid(grid* g, size_t n, size_t m, double v);
void copy_grid(grid* src, grid* dst);
void dump_grid(char* filename, grid *g);
void parse_file(char* filename, grid *g);
