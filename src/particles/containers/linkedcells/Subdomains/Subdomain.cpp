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
    //apply simple forces
    for (auto& cell : subdomainCells) {
        //first update position
        for (auto *particle: cell->getParticleReferences()) {
            const std::array<double, 3> new_v = particle->getV() + (delta_t / (2 * particle->getM())) * (particle->getF() + particle->getOldF());
            particle->setV(new_v);
        }
    }
    // apply pairwise forces



    // //TODO: port this to subdomains
    /*
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

    //#pragma omp parallel for schedule(dynamic)
    for (Cell* cell : occupied_cells_references) {
        // skip halo cells
        // if (cell->getCellType() == Cell::CellType::HALO) continue;

        for (auto it1 = cell->getParticleReferences().begin(); it1 != cell->getParticleReferences().end(); ++it1) {
            Particle* p = *it1;
            // calculate the forces between the particle and the particles in the same cell
            // uses direct sum with newtons third law
            for (auto it2 = (it1 + 1); it2 != cell->getParticleReferences().end(); ++it2) {
                Particle* q = *it2;
                std::array<double, 3> total_force{0, 0, 0};
                for (auto& force : force_sources) {
                    total_force = total_force + force->calculateForce(*p, *q);
                }
                p->setF(p->getF() + total_force);
                q->setF(q->getF() - total_force);
            }

            // calculate the forces between the particle and the particles in the neighbour cells
            for (Cell* neighbour : cell->getNeighbourReferences()) {
                if (cell->getAlreadyInfluencedBy().contains(neighbour)) continue;

                for (Particle* neighbour_particle : neighbour->getParticleReferences()) {
                    if (ArrayUtils::L2Norm(p->getX() - neighbour_particle->getX()) > cutoff_radius) continue;

                    for (const auto& force_source : force_sources) {
                        std::array<double, 3> force = force_source->calculateForce(*p, *neighbour_particle);
                        p->setF(p->getF() + force);
                        neighbour_particle->setF(neighbour_particle->getF() - force);
                    }
                }

                neighbour->addAlreadyInfluencedBy(cell); // this should be a loop carried dependency
            }
        }
    }

    // remove the periodic halo particles
    deleteHaloParticles();
    updateCellsParticleReferences();*/



    //upddate velocities
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
