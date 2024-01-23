#include "VerletFunctor.h"

#include "utils/ArrayUtils.h"

#include "omp.h"
#include "particles/containers/linkedcells/boundaries/ReflectiveBoundaryType.h"
#include "particles/containers/linkedcells/boundaries/OutflowBoundaryType.h"
#include "particles/containers/linkedcells/boundaries/PeriodicBoundaryType.h"

void VerletFunctor::step(std::unique_ptr<ParticleContainer> &particle_container,
                         const std::vector<std::shared_ptr<SimpleForceSource>> &simple_force_sources,
                         const std::vector<std::shared_ptr<PairwiseForceSource>> &pairwise_force_sources,
                         double delta_t) const {
    for (auto &p: *particle_container) {
        // update position
        const std::array<double, 3> new_x =
                p.getX() + delta_t * p.getV() + (delta_t * delta_t / (2 * p.getM())) * p.getF();
        p.setX(new_x);

        // reset forces
        p.setOldF(p.getF());
        p.setF({0, 0, 0});
    }

    /*
     * iteration over cells and give the cells
     */

    // calculate new forces
    particle_container->prepareForceCalculation();
    particle_container->applySimpleForces(simple_force_sources);
    particle_container->applyPairwiseForces(pairwise_force_sources);

    for (auto &p: *particle_container) {
        const std::array<double, 3> new_v = p.getV() + (delta_t / (2 * p.getM())) * (p.getF() + p.getOldF());
        p.setV(new_v);
    }

}

template<unsigned int N>
void VerletFunctor::templated_step(std::unique_ptr<ParticleContainer> &particle_container,
                                   const std::vector<std::shared_ptr<SimpleForceSource>> &simple_force_sources,
                                   const std::vector<std::shared_ptr<PairwiseForceSource>> &pairwise_force_sources,
                                   double delta_t) {
    static_assert(0 < N < 4);  // N should only be 1, 2 or 3

    if constexpr (N == 1) { // parallelization strategy 1: subdomains
        /*for (auto &subdomain : particle_container->) {

        }*/
        // TODO: maybe I can split it up into subdomains already here
        /*
         * for(subdomain: subdomains)
         *   subdomain.updatePositionSubdomain
         *   prepareForceCalculationSubdomain
         *
         *
         */
        auto numDomains = 10;
        std::map<unsigned int, Subdomain *> subdomainsStep = particle_container->getSubdomains();
        // update die positions
        // dann barrier
        // dann update den rest so gut es geht

#pragma omp parallel for schedule(static, 1)
        for (int i = 0; i < particle_container->getSubdomains().size(); ++i) {
            subdomainsStep.at(i)->updateParticlePositions();
        }
#pragma omp barrier
         particle_container->prepareForceCalculation();
         ReflectiveBoundaryType::applyBoundaryConditions(*particle_container);
         OutflowBoundaryType::applyBoundaryConditions(particle_container);
         PeriodicBoundaryType::applyBoundaryConditions(particle_container);

        // now the prepare force calculation steps
        /* particle_container->updatePositionSubdomain();
         particle_container->prepareForceCalculation(); // problem, here the positions will be updated, so the domains would have
         // to exchange information
         ReflectiveBoundaryType::applyBoundaryConditions(*this);
       OutflowBoundaryType::applyBoundaryConditions(*this);
    PeriodicBoundaryType::applyBoundaryConditions(*this);

    // clear the already influenced by vector in the cells
    // this is needed to prevent the two cells from affecting each other twice
    // since newtons third law is used
#pragma omp parallel for schedule(dynamic)
    for (Cell* cell : domain_cell_references) {
        cell->clearAlreadyInfluencedBy();
    }

         particle_container->applySimpleForcesDomains(simple_force_sources);
         particle_container->applyPairwiseForcesDomains(pairwise_force_sources);
         particle_container->updateVelocitySubdomain();*/

        particle_container->prepareForceCalculation();
#pragma omp parallel for schedule(static, 1)
        for (int i = 0; i < particle_container->getSubdomains().size(); ++i) {
            subdomainsStep.at(i)->updateSubdomain();
        }

    } else if (N == 2) { // parallelization strategy 2

    } else { // sequentiel implementaton

    }
}
