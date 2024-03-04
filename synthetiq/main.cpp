#include "algo.h"

/**
 * @brief The entry point of the program.
 * 
 * @param argc The number of command-line arguments.
 * @param argv An array of command-line arguments.
 * @return int The exit code of the program.
 */
int main(int argc, char* argv[]) {
    Algorithm algo = Algorithm(argc, argv);
    algo.run();
}
