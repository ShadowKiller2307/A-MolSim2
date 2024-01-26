#include "gtest/gtest.h"
#include "io/logger/Logger.h"
#include "particles/containers/linkedcells/LinkedCellsContainer.h"
#include "simulation/Simulation.h"
#include "simulation/SimulationUtils.h"
#include "utils/ArrayUtils.h"
#include "../../src/simulation/interceptors/diffusion/DiffusionInterceptor.h"


using BC = LinkedCellsContainer::BoundaryCondition;

/// @brief Test if the variance is calculated correctly for one particle after interacting with a periodic boundary
TEST(Diffusion,ProperDisplacementAfterPeriodic){
    Logger::logger->set_level(spdlog::level::err);

    std::vector<Particle> particles;
    std::array<double, 3> x1 = {5, 1, 1};
    std::array<double, 3> v1 = {1,0, 0};
    auto particle = Particle(x1, v1, 1, 0);

    particles.push_back(particle);
    particle.setDisplacementToAdd({0,0,0});
    SimulationParams params = TEST_DEFAULT_PARAMS_LENNARD_JONES;
    params.end_time = 1;
    params.delta_t = 0.1;

    params.container_type =
            SimulationParams::LinkedCellsType({6, 10, 10}, 3, {BC::PERIODIC, BC::PERIODIC,
                                                                 BC::OUTFLOW,BC::PERIODIC,
                                                                 BC::PERIODIC, BC::PERIODIC});
    Simulation simulation = Simulation(particles, params);
    DiffusionInterceptor interceptor = DiffusionInterceptor();
    interceptor.onSimulationStart(simulation);
    interceptor.setTestMode(true);
    auto overview = simulation.runSimulation();
    interceptor.operator()(1,simulation);

    EXPECT_NEAR(1.21,interceptor.getVariancesForTestMode()[0],1e-10);






}
