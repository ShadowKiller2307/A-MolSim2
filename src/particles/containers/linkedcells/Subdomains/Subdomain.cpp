#include "Subdomain.h"
#include <omp.h>
#include "utils/ArrayUtils.h"



Subdomain::Subdomain(double delta_t, double gravityConstant) : delta_t{delta_t}, gravityConstant{gravityConstant} {
    subdomainCells = {};
  /*  auto amountThreads = numThreads; //omp_get_num_threads();
    auto amountSubdomains = 10; // TODO: berechnen
    for (int i = 0; i < ; ++i) {

    }*/
}

//void Subdomain::initializeSubdomains() {}
void Subdomain::addCell(Cell* cellToAdd) {
    subdomainCells.emplace(cellToAdd);
}

void Subdomain::updateSubdomain() {
    /*for (auto& cell : subdomainCells) {
        //first update position
        for (auto* particle : cell->getParticleReferences()) {
            const std::array<double, 3> new_x =
                    particle->getX() + delta_t * particle->getV() + (delta_t * delta_t / (2 * particle->getM())) * particle->getF();
            particle->setX(new_x);
            // reset forces
            particle->setOldF(particle->getF());
            particle->setF({0, 0, 0});
        }
//#omp barrier // hier barrier
        // which thread updates the cells


        *//* prepare force calculation
        updateCellsParticleReferences();
        ReflectiveBoundaryType::pre(*this);
        OutflowBoundaryType::pre(*this);
        PeriodicBoundaryType::pre(*this);
        updateCellsParticleReferences();
         *//*

        // TODO: das hier für domains anpassen
       *//* particle_container->prepareForceCalculation();
        particle_container->applySimpleForces(simple_force_sources);
        particle_container->applyPairwiseForces(pairwise_force_sources);

        for (auto &p: *particle_container) {
            const std::array<double, 3> new_v = p.getV() + (delta_t / (2 * p.getM())) * (p.getF() + p.getOldF());
            p.setV(new_v);
        }
        *//*
    }*/


    for (auto& cell : subdomainCells) {
        //first update position
        for (auto *particle: cell->getParticleReferences()) {
            const std::array<double, 3> new_v = particle->getV() + (delta_t / (2 * particle->getM())) * (particle->getF() + particle->getOldF());
            particle->setV(new_v);
        }
    }
}

void Subdomain::updateParticlePositions() {
    for (auto& cell : subdomainCells) {
        //first update position
        for (auto *particle: cell->getParticleReferences()) {
            const std::array<double, 3> new_x =
                    particle->getX() + delta_t * particle->getV() +
                    (delta_t * delta_t / (2 * particle->getM())) * particle->getF();
            particle->setX(new_x);
            // reset forces
            particle->setOldF(particle->getF());
            particle->setF({0, 0, 0});
        }
    }
}
