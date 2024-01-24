#include "Subdomain.h"

Subdomain::Subdomain(double delta_t, double gravityConstant, double cutOffRadius) : delta_t{delta_t}, gravityConstant{gravityConstant}, cutoffRadius{cutOffRadius} {
    subdomainCells = {};
}

//void Subdomain::initializeSubdomains() {}
void Subdomain::addCell(Cell* cellToAdd) {
    subdomainCells.emplace(cellToAdd);
}

void Subdomain::updateSubdomain(const std::vector<std::shared_ptr<PairwiseForceSource>>& force_sources) {
    //apply simple forces
    for (auto& cell : subdomainCells) {
        for (auto* particle : cell->getParticleReferences()) { //directly apply gravitational
            //force here for performance
            auto currentF = particle->getF();
            currentF[1] += particle->getM() * (-gravityConstant);
            particle->setF(currentF);
        }
    }
    // apply pairwise forces
    //TODO: maybe domain occupied_cell_references
    for (Cell* cell : subdomainCells) {
        // skip halo cells
        // if (cell->getCellType() == Cell::CellType::HALO) continue;

        for (auto it1 = cell->getParticleReferences().begin(); it1 != cell->getParticleReferences().end(); ++it1) {
            Particle* p = *it1;
            // calculate the forces between the particle and the particles in the same cell
            // uses direct sum with newtons third law
            for (auto it2 = (it1 + 1); it2 != cell->getParticleReferences().end(); ++it2) {
                Particle* q = *it2;
                std::array<double, 3> total_force{0, 0, 0};
                for (auto& force : force_sources) {  //maybe inline function for the force ?
                    total_force = total_force + force->calculateForce(*p, *q);
                }
                p->setF(p->getF() + total_force);
                q->setF(q->getF() - total_force);
            }

            // calculate the forces between the particle and the particles in the neighbour cells
            //TODO: what to do if the neighbour cell is in a different subdomain?
            // -> is the already influenced by vector the best option?
            
            for (Cell* neighbour : cell->getNeighbourReferences()) {
                if (cell->getAlreadyInfluencedBy().contains(neighbour)) continue;

                /*for (Particle* neighbour_particle : neighbour->getParticleReferences()) {
                    if (ArrayUtils::L2Norm(p->getX() - neighbour_particle->getX()) > cutoffRadius) continue;

                    for (const auto& force_source : force_sources) {
                        std::array<double, 3> force = force_source->calculateForce(*p, *neighbour_particle);
                        p->setF(p->getF() + force);
                        neighbour_particle->setF(neighbour_particle->getF() - force);
                        // for cells at the border of subdomains maybe disable Newtons 3rd law
                        // or maybe lock the cells
                    }
                }

                neighbour->addAlreadyInfluencedBy(cell); // this should be a loop carried dependency
            }*/
        }
    }
    //update velocities
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
