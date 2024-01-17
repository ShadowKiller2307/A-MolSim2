#include "VerletFunctor.h"

#include "utils/ArrayUtils.h"

#include "omp.h"


void VerletFunctor::step(std::unique_ptr<ParticleContainer>& particle_container,
                         const std::vector<std::shared_ptr<SimpleForceSource>>& simple_force_sources,
                         const std::vector<std::shared_ptr<PairwiseForceSource>>& pairwise_force_sources, double delta_t) const {

   #pragma omp parallel for
    for (auto& p : *particle_container) {
        // update position
        const std::array<double, 3> new_x = p.getX() + delta_t * p.getV() + (delta_t * delta_t / (2 * p.getM())) * p.getF();
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

    // update velocity
    #pragma omp parallel for
    for (auto& p : *particle_container) {
        const std::array<double, 3> new_v = p.getV() + (delta_t / (2 * p.getM())) * (p.getF() + p.getOldF());
        p.setV(new_v);
    }
}
