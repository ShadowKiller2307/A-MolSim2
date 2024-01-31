#pragma once

#include <map>
#include <memory>
#include <vector>
#include <spdlog/fmt/chrono.h>

#include "integration/IntegrationMethods.h"
#include "particles/Particle.h"
#include "particles/containers/ParticleContainer.h"
#include "simulation/SimulationOverview.h"
#include "integration/methods/VerletFunctor.h"

class SimulationInterceptor;

class SimulationParams;

/**
 * @brief Class to run a simulation
 *
 * This class collects all the components needed to run a simulation, and provides a method to run it.
 */
/* template <unsigned N>*/
class Simulation {
public:
    /**
     * @brief Construct a new Simulation object and initialize all the necessary components
     *
     * @param particles Reference to the `ParticleContainer` on whose content the simulation is performed
     * @param params Parameters for the simulation. See the class `SimulationParams` for more information
     * @param integration_method The integration method to use for the simulation (Default: `IntegrationMethod::VERLET`)
     */
    Simulation(const std::vector<Particle> &particles, const SimulationParams &params,
               IntegrationMethod integration_method = IntegrationMethod::VERLET);

    ~Simulation();

    /**
     * @brief Runs the simulation, using the parameters given at construction and returns a `SimulationOverview` object containing some data
     *
     * @return SimulationOverview object containing some data about the simulation performed
     */
    //TODO: the templated function should be in the parent class, but what to do as
    // params will only be instantiated at runtime
    SimulationOverview runSimulation();

private:
/**
 * @brief Reference to the simulation parameters object
 */
const SimulationParams &params;

int strategy = 3;

//    static bool is_strategy_1;

/**
 * @brief Reference to the `ParticleContainer` on whose content the simulation is performed
 */
std::unique_ptr<ParticleContainer> particle_container;

std::unique_ptr<LinkedCellsContainer> linkedCellsContainer;

/**
 * @brief Functor used to integrate the particles
 */
std::unique_ptr<IntegrationFunctor> integration_functor;
// std::unique_ptr<VerletFunctor> verletFunctor;

std::vector<std::array<double, 3>> initial_pos_of_particles;

static void savePerformanceTest(const SimulationOverview &overview, const SimulationParams &params);

/**
 * Befriend the interceptors to allow them to access the private members of this class
 */
friend class ProgressBarInterceptor;

friend class FrameWriterInterceptor;

friend class ThermostatInterceptor;

friend class ParticleUpdateCounterInterceptor;

friend class RadialDistributionFunctionInterceptor;

friend class DiffusionInterceptor;

};