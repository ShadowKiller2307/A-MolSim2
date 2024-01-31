#include "ForcePicker.h"

#include <numeric>

#include "io/logger/Logger.h"
#include "physics/pairwiseforces/GravitationalForce.h"
#include "physics/pairwiseforces/LennardJonesForce.h"
#include "physics/pairwiseforces/HarmonicPotential.h"
#include "physics/simpleforces/GlobalDownwardsGravity.h"

const std::map<std::string, std::shared_ptr<SimpleForceSource>> get_supported_simple_forces()
{
    std::map<std::string, std::shared_ptr<SimpleForceSource>> force_names;

    auto global_downwards_gravity = std::make_shared<GlobalDownwardsGravity>(0);

    force_names.insert({std::string(*global_downwards_gravity), global_downwards_gravity});

    return force_names;
}

const std::map<std::string, std::shared_ptr<PairwiseForceSource>> get_supported_pairwise_forces()
{
    std::map<std::string, std::shared_ptr<PairwiseForceSource>> force_names;

    auto lennardjones = std::make_shared<LennardJonesForce>();
    auto gravitational = std::make_shared<GravitationalForce>();
    auto harmonic = std::make_shared<HarmonicPotential>(0.0, 0.0);

    force_names.insert({std::string(*lennardjones), lennardjones});
    force_names.insert({std::string(*gravitational), gravitational});
    force_names.insert({std::string(*harmonic), harmonic});

    return force_names;
}
