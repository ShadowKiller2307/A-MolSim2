#include "gtest/gtest.h"
#include "io/logger/Logger.h"
#include "particles/containers/linkedcells/LinkedCellsContainer.h"
#include "simulation/Simulation.h"
#include "simulation/SimulationUtils.h"
#include "utils/ArrayUtils.h"
#include "../../src/simulation/interceptors/diffusion/DiffusionInterceptor.h"
#include "io/input/xml/XMLFileReader.h"


using BC = LinkedCellsContainer::BoundaryCondition;

/// @brief Test if the variance is calculated correctly for one particle after interacting with a periodic boundary
TEST(Diffusion,ProperVarianceAfterPeriodic){
   Logger::logger->set_level(spdlog::level::err);


    XMLFileReader file_reader(true);

    auto [particles, params] = file_reader.readFile(FileLoader::get_input_file_path("test_diffusion.xml"));

   for(Particle& p:particles){

       p.setX({0.571774521296207, 2.356113654397094, 2.114357735420227});
   }

    auto conf = *params;


   Simulation simulation = Simulation(particles,conf);
   DiffusionInterceptor interceptor = DiffusionInterceptor();
   interceptor.onSimulationStart(simulation);
   std::vector<std::array<double,3>> previous_references = {{5,2,2}};
   interceptor.setTestMode(true);
   interceptor.setPreviousReferences(previous_references);
   interceptor.operator()(7,simulation);

   EXPECT_NEAR(interceptor.getVariancesForTestMode().at(0),2.610,1e-3);




}
