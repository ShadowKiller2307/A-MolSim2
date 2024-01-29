#include "gtest/gtest.h"
#include "io/logger/Logger.h"
#include "simulation/Simulation.h"
#include "simulation/SimulationUtils.h"
#include "utils/ArrayUtils.h"
#include "../../src/simulation/interceptors/radial_distribution_function/RadialDistributionFunctionInterceptor.h"
#include "io/input/xml/XMLFileReader.h"

/// @brief Tests if the bin_index and density are calculated correctly for one particle pair
TEST(RDF,ProperIndexAndDensity) {

    XMLFileReader file_reader(true);

    auto [particles, params] = file_reader.readFile(FileLoader::get_input_file_path("test_rdf.xml"));
    auto conf = *params;

    Simulation simulation = Simulation(particles,conf);
    RadialDistributionFunctionInterceptor interceptor = RadialDistributionFunctionInterceptor(1,20);
    interceptor.setTestMode(true);
    interceptor.onSimulationStart(simulation);

    interceptor.operator()(1,simulation);
    std::map<size_t, double> index_and_densities = interceptor.getDensitiesForTestMode();
    auto& pair = *index_and_densities.begin();


    EXPECT_EQ(pair.first,3);
    EXPECT_NEAR(pair.second,0.006452,1e-3);

}
