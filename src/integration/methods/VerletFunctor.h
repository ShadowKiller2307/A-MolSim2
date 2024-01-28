#pragma once

#include "integration/IntegrationFunctor.h"
#include "particles/containers/linkedcells/boundaries/ReflectiveBoundaryType.h"
#include "particles/containers/linkedcells/boundaries/OutflowBoundaryType.h"
#include "particles/containers/linkedcells/boundaries/PeriodicBoundaryType.h"


/**
 * @brief Implementation of the Strömer-Verlet integration method. Implements the interface IntegrationFunctor.
 *
 * Implements the IntegrationFunctor interface, and therefore updates all particles in the linkedCellsContainer according to the Verlet
 * integration method.
 */
class VerletFunctor : public IntegrationFunctor {
public:
    /**
     * @brief Performs one step with the respective integration method.
     *
     * @param particleContainer Container of particles on which the integration step is applied
     * @param simple_force_sources Vector of simple force sources which are used to calculate the new forces
     * @param pairwise_force_sources Vector of pairwise force sources which are used to calculate the new forces
     * @param delta_t Time step
     */
    void step(std::unique_ptr<ParticleContainer> &particleContainer,
              const std::vector<std::shared_ptr<SimpleForceSource>> &simple_force_sources,
              const std::vector<std::shared_ptr<PairwiseForceSource>> &pairwise_force_sources,
              double delta_t) const override;

    void parallel_step(
            std::unique_ptr<LinkedCellsContainer> &linkedCellsContainer, // maybe change this to LinkedCellsContainer
            const std::vector<std::shared_ptr<SimpleForceSource>> &simple_force_sources,
            const std::vector<std::shared_ptr<PairwiseForceSource>> &pairwise_force_sources, double delta_t,
            double gravityConstant, int strategy);
};