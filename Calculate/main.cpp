
#include "Solver.h"


int main(int argc, char *argv[])
{
    Solve::PoisonEquation("../Mesh/Drop.msh", argc, argv);
    return 0;
}