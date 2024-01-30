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
//    SimulationOverview runSimulation() {
//        size_t iteration = 0;
//        double simulated_time = 0;
//
//        // Calculate initial forces
//        particle_container->prepareForceCalculation();
//        particle_container->applySimpleForces(params.simple_forces);
//        particle_container->applyPairwiseForces(params.pairwise_forces);
//
//        Logger::logger->info("Simulation started...");
//
//        std::time_t t_start_helper = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
//        Logger::logger->info("Start time: {}", fmt::format("{:%A %Y-%m-%d %H:%M:%S}", fmt::localtime(t_start_helper)));
//
//        // Notify interceptors that the simulation is about to start
//        for (auto &interceptor: params.interceptors) {
//            (*interceptor).onSimulationStart(*this);
//        }
//
//        auto t_start = std::chrono::high_resolution_clock::now();
//
//        *//*if (N == 1 || N == 2) {
//                static_cast<LinkedCellsContainer>(*particle_container);
//        // }*//*
//
//        std::unique_ptr<VerletFunctor> verletFunctor;
//
//        while (simulated_time < params.end_time) {
//            verletFunctor->templated_step<N>();
//            *//* if constexpr (N == 1) { //
//                    verletFunctor->templated_step<N>();
//            //integration_functor->step(particle_container, params.simple_forces, params.pairwise_forces, params.delta_t);
//        }
//        else if constexpr (N == 2) {
//
//        } else {
//
//        }
//        *//*
//                ++iteration;
//        simulated_time += params.delta_t;
//
//        // Notify interceptors of the current iteration
//        for (auto &interceptor: params.interceptors) {
//            (*interceptor).notify(iteration, *this);
//        }
//    }
//
//    auto t_end = std::chrono::high_resolution_clock::now();
//
//    // Notify interceptors that the simulation has ended
//    for (
//    auto &interceptor: params.interceptors
//    ) {
//        (*interceptor).onSimulationEnd(iteration, *this);
//    }
//
//    Logger::logger->info("Simulation finished.");
//    Logger::logger->info("End time: {}", fmt::format("{:%A %Y-%m-%d %H:%M:%S}",
//    fmt::localtime(t_end)
//    ));
//
//    std::vector<std::string> interceptor_summaries;
//    for (
//    auto &interceptor: params.interceptors
//    ) {
//        auto summary = std::string(*interceptor);
//        if (!summary.empty()) {
//            interceptor_summaries.push_back(summary);
//        }
//    }
//
//    auto total_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();
//
//    SimulationOverview overview{params, total_time_ms / 1000.0, iteration, interceptor_summaries,
//                                std::vector<Particle>(particle_container->begin(), particle_container->end())};
//
//    if (params.performance_test) {
//        savePerformanceTest(overview, params);
//    }
//
//    return
//    overview;
//
//}

//    static void setStrategy1(bool to_use);
//    static bool getIsStrategy1();

private:
/**
 * @brief Reference to the simulation parameters object
 */
const SimulationParams &params;

int strategy = 1;

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