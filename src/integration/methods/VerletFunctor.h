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
    void step(std::unique_ptr<ParticleContainer>& particleContainer,
              const std::vector<std::shared_ptr<SimpleForceSource>>& simple_force_sources,
              const std::vector<std::shared_ptr<PairwiseForceSource>>& pairwise_force_sources, double delta_t) const override;

    template<unsigned N>
    void templated_step(std::unique_ptr<LinkedCellsContainer>& linkedCellsContainer, // maybe change this to LinkedCellsContainer
                                 const std::vector<std::shared_ptr<SimpleForceSource>>& simple_force_sources,
                                 const std::vector<std::shared_ptr<PairwiseForceSource>>& pairwise_force_sources, double delta_t) {
      /*  static_assert(0 < N < 4);  // N should only be 1, 2 or 3

        if constexpr (N == 1) { // parallelization strategy 1: subdomains
            *//*for (auto &subdomain : linkedCellsContainer->) {

            }*//*
            auto numDomains = 10;
            std::map<unsigned int, Subdomain *> subdomainsStep = linkedCellsContainer->getSubdomains();
            // update die positions
            // dann barrier
            // dann update den rest so gut es geht

#pragma omp parallel for schedule(static, 1)
//TODO: maybe tasks are better here
            for (int i = 0; i < linkedCellsContainer->getSubdomains().size(); ++i) {
                subdomainsStep.at(i)->updateParticlePositions();
            }
#pragma omp barrier
#pragma omp parallel for schedule(static, 1)
            linkedCellsContainer->prepareForceCalculation();

            ReflectiveBoundaryType::applyBoundaryConditions(*linkedCellsContainer);
            OutflowBoundaryType::applyBoundaryConditions(*linkedCellsContainer);
            PeriodicBoundaryType::applyBoundaryConditions(*linkedCellsContainer);
#pragma omp parallel for schedule(static, 1)
#pragma omp single
            //TODO: only one thread apply the boundary conditions

            for (int i = 0; i < linkedCellsContainer->getSubdomains().size(); ++i) {
                subdomainsStep.at(i)->updateSubdomain(pairwise_force_sources);
            }

        } else if (N == 2) { // parallelization strategy 2

        } else { // sequentiel implementaton, maybe call the sequential implementation only
            // if the old step method is called and not the templated step

        }*/
    }

/*private:
    unsigned parallelization strategy*/
};