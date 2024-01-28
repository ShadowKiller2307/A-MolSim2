#include "io/cli/CLIParser.h"
#include "io/input/FileInputHandler.h"
#include "simulation/Simulation.h"
#include "simulation/SimulationParams.h"

int main(int argc, char* argsv[]) {
    // Parse CLI arguments
    CLIParams params_cli = parse_arguments(argc, argsv);

    // Parse input file
    auto [initial_particles, simulation_arguments] = FileInputHandler::readFile(params_cli.input_file_path, params_cli.fresh);

    // Combine parameters from CLI and input file
    SimulationParams params = merge_parameters(params_cli, simulation_arguments);
    params.num_particles = initial_particles.size();

#ifdef OMP_NUM_THREADS
    omp_set_num_threads(OMP_NUM_THREADS); // set the number of threads to the amount the user specifies
#endif

#ifdef PARTICLES
    std::cout << "ifdef particles funktioniert" << std::endl;
#endif

#ifdef SUBDOMAIN
    std::cout << "ifdef subdomain funktioniert" << std::endl;
#endif
    // Initialize simulation
    Simulation simulation{initial_particles, params};

    // Print simulation info
    params.logSummary();

    auto overview = simulation.runSimulation(); //<3>();

    // Print simulation overview
    overview.logSummary();

    return 0;
}
