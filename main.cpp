#include "utils.hpp"
#include "GeneticAlgorithm.hpp"

int main(int argc, char* argv[])
{
    Config cf;

    if (!parsing(argc, argv, cf))
        return 0;

    GeneticAlgorithm* GA = new GeneticAlgorithm();
    GA->initialize(cf);
    GA->run();
    GA->output_ith_best_layout();

    return 0;
}