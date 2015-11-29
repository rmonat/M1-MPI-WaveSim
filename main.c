// TODO : iteration et steps plutôt des unsigned int
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "automaton.h"
#include "args.h"
#include "step0.h"


int main(int argc, char** argv)
{
    inst instance;
    parse_options(&instance, argc, argv);

    if(instance.step == 0)
    {
	grid g;
	parse_file(instance.input_path, &g);
	step0(instance, &g);
    }
    return 0;
}

