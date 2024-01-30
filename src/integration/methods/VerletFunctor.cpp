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

    // calculate new forces
    //std::cout << "begin calculate new forces" << std::endl;
    particle_container->prepareForceCalculation();  // here the boundary condition pre functions will be called
    particle_container->applySimpleForces(simple_force_sources);
    particle_container->applyPairwiseForces(pairwise_force_sources); // here the boundary condition apply
    //   std::cout << "end calculate new forces" << std::endl;
    // boundary functions will be called

    // std::cout << "begin velocity calculation begin: " << std::endl;
    for (auto &p: *particle_container) {
        const std::array<double, 3> new_v = p.getV() + (delta_t / (2 * p.getM())) * (p.getF() + p.getOldF());
        p.setV(new_v);
    }
    // std::cout << "verlet functor step end" << std::endl;

}

void VerletFunctor::parallel_step(std::unique_ptr<ParticleContainer> &linkedCellsContainer,
                                  const std::vector<std::shared_ptr<SimpleForceSource>> &simple_force_sources,
                                  const std::vector<std::shared_ptr<PairwiseForceSource>> &pairwise_force_sources,
                                  double delta_t, double gravityConstant, int strategy) {
   /* if (strategy == 1) { // parallelization strategy 1: subdomains
        auto numDomains = 10;
        std::vector<Subdomain *> subdomainsStep = linkedCellsContainer->getSubdomainsVector();
        // update die positions
        // dann barrier
        // dann update den rest so gut es geht
#pragma omp parallel
        {
            auto numThreads = omp_get_num_threads();
            std::cout << "Number of threads" << std::endl;
#pragma omp for schedule(static, 1) // static, 1 because every thread should work on a subdomain
            for (int i = 0; i < linkedCellsContainer->getSubdomains().size(); ++i) {
                subdomainsStep.at(i)->updateParticlePositions();
            }  // <-- omp puts barrier automatically here
            linkedCellsContainer->prepareForceCalculation();
            ReflectiveBoundaryType::applyBoundaryConditions(*linkedCellsContainer);
            OutflowBoundaryType::applyBoundaryConditions(*linkedCellsContainer);
            PeriodicBoundaryType::applyBoundaryConditions(*linkedCellsContainer);

            /// force calculations
#pragma omp for schedule(static, 1)
            for (int i = 0; i < linkedCellsContainer->getSubdomains().size(); ++i) {
                subdomainsStep.at(i)->updateSubdomain(pairwise_force_sources);
            }
        }

        linkedCellsContainer->deleteHaloParticles();
        linkedCellsContainer->updateCellsParticleReferences();

    } else if (strategy == 2) { // parallelization strategy 2 : particles parallelization
        *//**
         * @brief: for parallelizaton strategy 2, each cell will be processed sequentielly,
         * within the cells threads are spawned to update the particles
         *//*
        //std::cout << "reach parallel step" << std::endl;
*//*#pragma omp parallel
        {*//*
        /// update the positions
        //   std::cout << "Strategy particles" << std::endl;
#pragma omp parallel for schedule(dynamic)
        for (auto &p: *linkedCellsContainer) {
            // std::cout << omp_get_num_threads() << std::endl; // at the moment only one thread
            // update position
            const std::array<double, 3> new_x =
                    p.getX() + delta_t * p.getV() + (delta_t * delta_t / (2 * p.getM())) * p.getF();
            p.setX(new_x);

            // reset forces
            p.setOldF(p.getF());
            p.setF({0, 0, 0});
        }
   //     std::cout << "after position calculation" << std::endl;

#pragma omp single
        {
            linkedCellsContainer->prepareForceCalculation();
        }

        /// apply simple forces
#pragma omp parallel for schedule(dynamic)
        for (auto &p: *linkedCellsContainer) { //directly apply gravitational
            //force here for performance
            auto currentF = p.getF();
            currentF[1] += p.getM() * (-gravityConstant);
            p.setF(currentF);
        }
       // std::cout << "after simple force calculation" << std::endl;

        /// apply pairwise forces

*//*#pragma omp single
            {*//*
        linkedCellsContainer->applyPairwiseForcesOptimized(pairwise_force_sources);
//        std::cout << "after pairwise force calculation" << std::endl;
        // in deleteHalo or updateParticleReferences a double free was cause by two threads simultanously trying to
        // free the same particle pointer
        //     }

        //update the velocity
#pragma omp parallel for schedule(dynamic)
        for (auto &p: *linkedCellsContainer) {
            const std::array<double, 3> new_v = p.getV() + (delta_t / (2 * p.getM())) * (p.getF() + p.getOldF());
            p.setV(new_v);
        }
        //}
    }

    else if (strategy == 3) {
#pragma omp parallel
        {
            // parallelisation over particles without particle locks


        }
    } else {
        return;
    }*/
}