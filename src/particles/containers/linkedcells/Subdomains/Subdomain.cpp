#include "Subdomain.h"
#include "particles/containers/linkedcells/LinkedCellsContainer.h"


Subdomain::Subdomain(double delta_t, double gravityConstant, double cutOffRadius,
                     std::unique_ptr<LinkedCellsContainer> linkedCellsContainer
) : delta_t{delta_t}, gravityConstant{gravityConstant}, cutoffRadius{cutOffRadius},
    linkedCellsContainer{std::move(linkedCellsContainer)} {
    subdomainCells = {};
}

//void Subdomain::initializeSubdomains() {}
void Subdomain::addCell(bool isAtSubdomainBorder, Cell *cellToAdd) {
    subdomainCells.emplace(isAtSubdomainBorder, cellToAdd);
}

void Subdomain::calculateForcesBetweenCells(Cell *one, Cell *two) {

}

//TODO: need to change the order of acquiring the locks because
// at the moment deadlocks can occur
void Subdomain::updateSubdomain(const std::vector<std::shared_ptr<PairwiseForceSource>> &force_sources) {
    //apply simple forces
    for (auto &cell: subdomainCells) {
        for (auto *particle: cell.second->getParticleReferences()) { //directly apply gravitational
            //force here for performance
            auto currentF = particle->getF();
            currentF[1] += particle->getM() * (-gravityConstant);
            particle->setF(currentF);
        }
    }
    // apply pairwise forces
    //TODO: maybe domain occupied_cell_references
    for (auto &cell: subdomainCells) {
        /**
        only lock the cells which are at domain borders, because these are the only cells
        where race conflicts can occur(because every domain has at most one thread, so only one thread will
        work on the inner subdomain cells)
        whether the cell is at the subdomain border can be deduced from cell.first
        //
         */
        if (cell.first) {
            omp_set_lock(cell.second->getLock());
        }
        for (auto it1 = cell.second->getParticleReferences().begin();
             it1 != cell.second->getParticleReferences().end(); ++it1) {
            Particle *p = *it1;
            // calculate the forces between the particle and the particles in the same cell
            // uses direct sum with newtons third law
            for (auto it2 = (it1 + 1); it2 != cell.second->getParticleReferences().end(); ++it2) {
                Particle *q = *it2;
                std::array<double, 3> total_force{0, 0, 0};
                for (auto &force: force_sources) {  //maybe inline function for the force ?
                    total_force = total_force + force->calculateForce(*p, *q);
                }
                p->setF(p->getF() + total_force);
                q->setF(q->getF() - total_force);
            }
        }
        if (cell.first) {
            omp_set_lock(cell.second->getLock());
        }
        // calculate the forces between the particles in the current cell
        // and particles in the neighbouring cells
        for (Cell *neighbour: cell.second->getNeighboursToComputeForcesWith()) {
            //if (neighbour.first) {
            /// think a global lock order is established because of the deterministic force calculation
            omp_set_lock(cell.second->getLock());
            omp_set_lock(neighbour->getLock());
            //}
            for (Particle *p: cell.second->getParticleReferences()) {
                for (Particle *neighbour_particle: neighbour->getParticleReferences()) {
                    if (ArrayUtils::L2Norm(p->getX() - neighbour_particle->getX()) > cutoffRadius) continue;
                    for (const auto &force_source: force_sources) {
                        std::array<double, 3> force = force_source->calculateForce(*p, *neighbour_particle);
                        p->setF(p->getF() + force);
                        neighbour_particle->setF(neighbour_particle->getF() - force);
                    }
                }
            }
            // free the lock for the neighbour
            // if (cell.first) {
            omp_unset_lock(cell.second->getLock());
            omp_unset_lock(neighbour->getLock());
            //}
        }

        /*  // free the lock for the current cell
          if (cell.first) {
              omp_unset_lock(cell.second->getLock());
          }*/

    }
//update velocities
    for (
        auto &cell
            : subdomainCells) {
        for (
            auto *particle
                : cell.second->
                getParticleReferences()
                ) {
            const std::array<double, 3> new_v =
                    particle->getV() + (delta_t / (2 * particle->getM())) * (particle->getF() + particle->getOldF());
            particle->
                    setV(new_v);
        }
    }
}

void Subdomain::updateParticlePositions() {
    for (auto &cell: subdomainCells) {
        //first update position
        for (auto *particle: cell.second->getParticleReferences()) {
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
