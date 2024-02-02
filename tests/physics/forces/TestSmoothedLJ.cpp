#include <gtest/gtest.h>

#include "physics/pairwiseforces/SmoothedLJ.h"
#include "utils/ArrayUtils.h"
#include "io/logger/Logger.h"

///@brief Tests the smoothed LJ calculation for the middle case of the formula
TEST(SmoothedLJ,ProperCalculationMiddleCase){

    auto p1 = Particle({3,5,6.2},{0,0,0},1,0,1,1);
    auto p2 = Particle({4.76,4.16,5.43},{0,0,0},1,1,1,1);

    SmoothedLJ smoothed_lj = SmoothedLJ(1.9,2.3);
    std::array<double,3> result = smoothed_lj.calculateForce(p1,p2);
    std::cout<<"The result: "<< result << std::endl;
    Logger::logger->error("The result: ", result[0]);
    EXPECT_NEAR(result[0],0.000831,1e-3);
    EXPECT_NEAR(result[1],-0.000396,1e-3);
    EXPECT_NEAR(result[2],-0.000363,1e-3);

}



