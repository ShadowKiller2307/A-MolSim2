#include "LinkedCellsContainer.h"

#include <cmath>

#include "cells/Cell.h"
#include "io/logger/Logger.h"
#include "particles/containers/linkedcells/boundaries/OutflowBoundaryType.h"
#include "particles/containers/linkedcells/boundaries/PeriodicBoundaryType.h"
#include "particles/containers/linkedcells/boundaries/ReflectiveBoundaryType.h"
#include "physics/pairwiseforces/LennardJonesForce.h"
#include "utils/ArrayUtils.h"
#include <omp.h>

/*
    Methods of the LinkedCellsContainer
*/
LinkedCellsContainer::LinkedCellsContainer(const std::array<double, 3> &_domain_size, double _cutoff_radius,
                                           const std::vector<std::shared_ptr<PairwiseForceSource>> &pairwiseForceSources,
                                           const std::array<BoundaryCondition, 6> &_boundary_types, int _n,
                                           int numThreads)
        : domain_size(_domain_size), cutoff_radius(_cutoff_radius), pairwise_force_sources(pairwiseForceSources),
          boundary_types(_boundary_types), numThreads(numThreads) {

    // calculate the number of cells in each dimension
    domain_num_cells = {std::max(static_cast<int>(std::floor(_domain_size[0] / cutoff_radius)), 1),
                        std::max(static_cast<int>(std::floor(_domain_size[1] / cutoff_radius)), 1),
                        std::max(static_cast<int>(std::floor(_domain_size[2] / cutoff_radius)), 1)};

    cell_size = {_domain_size[0] / domain_num_cells[0], _domain_size[1] / domain_num_cells[1],
                 _domain_size[2] / domain_num_cells[2]};

    // reserve the memory for the cells
    cells.reserve((domain_num_cells[0] + 2) * (domain_num_cells[1] + 2) * (domain_num_cells[2] + 2));

    // create the cells with the correct cell-type and add them to the cells vector and the corresponding cell reference vector
    initCells();

    // add the neighbour references to the cells
    initCellNeighbourReferences();

    int numThreadsInit;
    if(numThreads == 8) {
        numThreadsInit = numThreads*2;
    }
    else {
        numThreadsInit = numThreads;
    }
    initSubdomains(numThreadsInit); // initialize the subdomains for parallelization strategy 2

    // reserve the memory for the particles to prevent reallocation during insertion
    particles.reserve(_n);

    Logger::logger->debug("Created LinkedCellsContainer with boundaries [{}, {}, {}] and cutoff radius {}",
                          domain_size[0], domain_size[1],
                          domain_size[2], cutoff_radius);
    Logger::logger->debug(
            "Created LinkedCellsContainer with {} domain cells (of which {} are at the boundary) and {} halo cells",
            domain_cell_references.size(), boundary_cell_references.size(), halo_cell_references.size());
    Logger::logger->debug("Cells per dimension: [{}, {}, {}]", domain_num_cells[0], domain_num_cells[1],
                          domain_num_cells[2]);
    Logger::logger->debug("Calculated cell size: [{}, {}, {}]", cell_size[0], cell_size[1], cell_size[2]);
}

void LinkedCellsContainer::addParticle(const Particle &p) {
    Cell *cell = particlePosToCell(p.getX());

    if (cell == nullptr) {
        Logger::logger->error("Particle to insert is out of bounds, position: [{}, {}, {}]", p.getX()[0], p.getX()[1],
                              p.getX()[2]);
        throw std::runtime_error("Attempted to insert particle out of bounds");
    }
    size_t old_capacity = particles.capacity();
    particles.push_back(p);

    if (old_capacity != particles.capacity()) {
        updateCellsParticleReferences();
    } else {
        cell->addParticleReference(&particles.back());
    }
}

void LinkedCellsContainer::addParticle(Particle &&p) {
    Cell *cell = particlePosToCell(p.getX());

    if (cell == nullptr) {
        Logger::logger->error("Particle to insert is outside of cells. Position: [{}, {}, {}]", p.getX()[0],
                              p.getX()[1], p.getX()[2]);
        throw std::runtime_error("Attempted to insert particle out of bounds");
    }

    size_t old_capacity = particles.capacity();
    particles.push_back(std::move(p));

    if (old_capacity != particles.capacity()) {
        updateCellsParticleReferences();
    } else {
        cell->addParticleReference(&particles.back());
    }
}

void LinkedCellsContainer::prepareForceCalculation() {
    // update the particle references in the cells in case the particles have moved
    updateCellsParticleReferences();

    //ReflectiveBoundaryType::pre(*this); //does nothing
    OutflowBoundaryType::pre(*this);   // deletes the halo cells
    PeriodicBoundaryType::pre(*this);  // moves the particles which left the domain to the other side

    // update the particle references in the cells in case the particles have moved
    updateCellsParticleReferences();
}

void
LinkedCellsContainer::applySimpleForces(const std::vector<std::shared_ptr<SimpleForceSource>> &simple_force_sources) {
    // std::cout << "apply simple forces begin" << std::endl;
    // strategy 1: distribute the particles over the threads
    //#pragma omp parallel for schedule(static, 100)
    //#pragma omp parallel for schedule(dynamic)
    for (Particle &p: particles) {
        for (const auto &force_source: simple_force_sources) {
            p.setF(p.getF() + force_source->calculateForce(p));
        }
    }
}

//TODO: there is a deadlock inside this function
void LinkedCellsContainer::applyPairwiseForcesOptimized(
        const std::vector<std::shared_ptr<PairwiseForceSource>> &force_sources) {
#ifdef _OPENMP
    //  std::cout << "reach begin of apply pairwise forces optimized" << std::endl;
/*#pragma omp single
    {*/
    ReflectiveBoundaryType::applyBoundaryConditions(*this);
    OutflowBoundaryType::applyBoundaryConditions(*this);
    PeriodicBoundaryType::applyBoundaryConditions(*this);

    for (Cell *cell: domain_cell_references) {
        cell->clearAlreadyInfluencedBy();
    }

    for (Cell *cell: occupied_cells_references) {
        // skip halo cells
        if (cell->getCellType() == Cell::CellType::HALO) continue;
        size_t amountOfParticles = cell->getParticleReferences().size();
#pragma omp parallel
        {
#pragma omp for schedule(dynamic)
            /// at the moment, the particles get locked
            for (size_t i = 0; i < amountOfParticles; ++i) {
                if (i >= cell->getParticleReferences().size()) {
                    throw std::out_of_range("Wrong particle access");
                }

                Particle *p = cell->getParticleReferences().at(i);
                for (size_t j = i + 1; j < amountOfParticles; ++j) {
                    //  std::cout << " inner force calculation begin" << std::endl;
                    if (j >= cell->getParticleReferences().size()) {
                        throw std::out_of_range("Wrong particle access");
                    }
                    Particle *q = cell->getParticleReferences().at(j);
                    std::array<double, 3> total_force{0, 0, 0};
                    for (auto &force: force_sources) {
                        total_force = total_force + force->calculateForce(*p, *q);
                    }
                    omp_set_lock(p->getLock());
                    p->setF(p->getF() + total_force);
                    omp_unset_lock(p->getLock());
                    omp_set_lock(q->getLock());
                    q->setF(q->getF() - total_force);
                    omp_unset_lock(q->getLock());
                    //std::cout << " inner force calculation end" << std::endl;
                }
                // std::cout << "before cell neighbours" << std::endl;
                for (Cell *neighbour: cell->getNeighboursToComputeForcesWith()) {
                    // std::cout << " outer force calculation begin" << std::endl;
                    //   if (cell->getAlreadyInfluencedBy().contains(neighbour)) continue;
                    for (Particle *neighbour_particle: neighbour->getParticleReferences()) {
                        if (ArrayUtils::L2Norm(p->getX() - neighbour_particle->getX()) > cutoff_radius) continue;

                        for (const auto &force_source: force_sources) {
                            std::array<double, 3> force = force_source->calculateForce(*p, *neighbour_particle);
                            omp_set_lock(p->getLock());
                            p->setF(p->getF() + force);
                            omp_unset_lock(p->getLock());
                            omp_set_lock(neighbour_particle->getLock());
                            neighbour_particle->setF(neighbour_particle->getF() - force);
                            omp_unset_lock(neighbour_particle->getLock());
                        }
                    }
                }

            }
        }
    }
    deleteHaloParticles();
    updateCellsParticleReferences();
#endif
}


void LinkedCellsContainer::applyPairwiseForces(const std::vector<std::shared_ptr<PairwiseForceSource>> &force_sources) {
    ReflectiveBoundaryType::applyBoundaryConditions(*this);
    OutflowBoundaryType::applyBoundaryConditions(*this);
    PeriodicBoundaryType::applyBoundaryConditions(*this);

    // clear the already influenced by vector in the cells
    // this is needed to prevent the two cells from affecting each other twice
    // since newtons third law is used
    for (Cell *cell: domain_cell_references) {
        cell->clearAlreadyInfluencedBy();
    }

    for (Cell *cell: occupied_cells_references) {
        // skip halo cells
        if (cell->getCellType() == Cell::CellType::HALO) continue;

        for (auto it1 = cell->getParticleReferences().begin(); it1 != cell->getParticleReferences().end(); ++it1) {
            Particle *p = *it1;
            // calculate the forces between the particle and the particles in the same cell
            // uses direct sum with newtons third law
            for (auto it2 = (it1 + 1); it2 != cell->getParticleReferences().end(); ++it2) {
                Particle *q = *it2;
                std::array<double, 3> total_force{0, 0, 0};
                for (auto &force: force_sources) {
                    total_force = total_force + force->calculateForce(*p, *q);
                }
                p->setF(p->getF() + total_force);
                q->setF(q->getF() - total_force);
            }

            // calculate the forces between the particle and the particles in the neighbour cells
            for (Cell *neighbour: cell->getNeighbourReferences()) {
                if (cell->getAlreadyInfluencedBy().contains(neighbour)) continue;

                for (Particle *neighbour_particle: neighbour->getParticleReferences()) {
                    if (ArrayUtils::L2Norm(p->getX() - neighbour_particle->getX()) > cutoff_radius) continue;

                    for (const auto &force_source: force_sources) {
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
    updateCellsParticleReferences();
}

void LinkedCellsContainer::reserve(size_t n) {
    Logger::logger->debug("Reserving space for {} particles", n);

    size_t old_capacity = particles.capacity();
    particles.reserve(n);

    if (old_capacity != particles.capacity()) {
        updateCellsParticleReferences();
    }
}

size_t LinkedCellsContainer::size() const { return particles.size(); }

Particle &LinkedCellsContainer::operator[](int i) { return particles[i]; }

std::vector<Particle>::iterator LinkedCellsContainer::begin() { return particles.begin(); }

std::vector<Particle>::iterator LinkedCellsContainer::end() { return particles.end(); }

std::vector<Particle>::const_iterator LinkedCellsContainer::begin() const { return particles.begin(); }

std::vector<Particle>::const_iterator LinkedCellsContainer::end() const { return particles.end(); }

const std::vector<Particle> &LinkedCellsContainer::getParticles() const { return particles; }

const std::array<double, 3> &LinkedCellsContainer::getDomainSize() const { return domain_size; }

double LinkedCellsContainer::getCutoffRadius() const { return cutoff_radius; }

const std::vector<Cell> &LinkedCellsContainer::getCells() { return cells; }

const std::vector<Cell *> &LinkedCellsContainer::getBoundaryCells() const { return boundary_cell_references; }

const std::array<double, 3> &LinkedCellsContainer::getCellSize() const { return cell_size; }

const std::array<int, 3> &LinkedCellsContainer::getDomainNumCells() const { return domain_num_cells; }

int LinkedCellsContainer::cellCoordToCellIndex(int cx, int cy, int cz) const {
    if (cx < -1 || cx > domain_num_cells[0] || cy < -1 || cy > domain_num_cells[1] || cz < -1 ||
        cz > domain_num_cells[2]) {
        return -1;
    }
    return (cx + 1) * (domain_num_cells[1] + 2) * (domain_num_cells[2] + 2) + (cy + 1) * (domain_num_cells[2] + 2) +
           (cz + 1);
}

int LinkedCellsContainer::findCellForParticle(const std::array<double, 3> &pos) {
    int cx = static_cast<int>(std::floor(pos[0] / cell_size[0]));
    int cy = static_cast<int>(std::floor(pos[1] / cell_size[1]));
    int cz = static_cast<int>(std::floor(pos[2] / cell_size[2]));

    return cellCoordToCellIndex(cx, cy, cz);
}

Cell *LinkedCellsContainer::particlePosToCell(const std::array<double, 3> &pos) {
    return particlePosToCell(pos[0], pos[1], pos[2]);
}

Cell *LinkedCellsContainer::particlePosToCell(double x, double y, double z) {
    int cx = static_cast<int>(std::floor(x / cell_size[0]));
    int cy = static_cast<int>(std::floor(y / cell_size[1]));
    int cz = static_cast<int>(std::floor(z / cell_size[2]));

    int cell_index = cellCoordToCellIndex(cx, cy, cz);
    if (cell_index == -1) {
        Logger::logger->error("Particle is outside of cells. Position: [{}, {}, {}]", x, y, z);
        throw std::runtime_error("A particle is outside of the cells");
    }

    return &cells[cell_index];
}

std::string LinkedCellsContainer::boundaryConditionToString(const BoundaryCondition &bc) {
    switch (bc) {
        case BoundaryCondition::OUTFLOW:
            return "Outflow";
        case BoundaryCondition::REFLECTIVE:
            return "Reflective";
        case BoundaryCondition::PERIODIC:
            return "Periodic";
        default:
            return "Unknown";
    }
}

/*
    Private methods of the LinkedCellsContainer
*/

void LinkedCellsContainer::initCells() {
    for (int cx = -1; cx < domain_num_cells[0] + 1; ++cx) {
        for (int cy = -1; cy < domain_num_cells[1] + 1; ++cy) {
            for (int cz = -1; cz < domain_num_cells[2] + 1; ++cz) {
                if (cx < 0 || cx >= domain_num_cells[0] || cy < 0 || cy >= domain_num_cells[1] || cz < 0 ||
                    cz >= domain_num_cells[2]) {
                    Cell new_cell(Cell::CellType::HALO, cells.size());
                    cells.push_back(new_cell);
                    halo_cell_references.push_back(&cells.back());

                    if (cx == -1) {
                        left_halo_cell_references.push_back(&cells.back());
                    }
                    if (cx == domain_num_cells[0]) {
                        right_halo_cell_references.push_back(&cells.back());
                    }
                    if (cy == -1) {
                        bottom_halo_cell_references.push_back(&cells.back());
                    }
                    if (cy == domain_num_cells[1]) {
                        top_halo_cell_references.push_back(&cells.back());
                    }
                    if (cz == -1) {
                        back_halo_cell_references.push_back(&cells.back());
                    }
                    if (cz == domain_num_cells[2]) {
                        front_halo_cell_references.push_back(&cells.back());
                    }
                } else if (cx == 0 || cx == domain_num_cells[0] - 1 || cy == 0 || cy == domain_num_cells[1] - 1 ||
                           cz == 0 ||
                           cz == domain_num_cells[2] - 1) {
                    Cell new_cell(Cell::CellType::BOUNDARY, cells.size());
                    cells.push_back(new_cell);
                    boundary_cell_references.push_back(&cells.back());
                    domain_cell_references.push_back(&cells.back());

                    if (cx == 0) {
                        left_boundary_cell_references.push_back(&cells.back());
                    }
                    if (cx == domain_num_cells[0] - 1) {
                        right_boundary_cell_references.push_back(&cells.back());
                    }
                    if (cy == 0) {
                        bottom_boundary_cell_references.push_back(&cells.back());
                    }
                    if (cy == domain_num_cells[1] - 1) {
                        top_boundary_cell_references.push_back(&cells.back());
                    }
                    if (cz == 0) {
                        back_boundary_cell_references.push_back(&cells.back());
                    }
                    if (cz == domain_num_cells[2] - 1) {
                        front_boundary_cell_references.push_back(&cells.back());
                    }
                } else {
                    Cell new_cell(Cell::CellType::INNER, cells.size());
                    cells.push_back(new_cell);
                    domain_cell_references.push_back(&cells.back());
                }
            }
        }
    }
}

std::array<double, 3> LinkedCellsContainer::computeSubdomainsPerDimension(int numThreads) {
    // the domain can be divided into 1, 2, 4, 8, 14, 16, 28 and 56 subdomains
    std::array<int, 8> numThreadOptions{1, 2, 4, 8, 14, 16, 28, 56};
    auto numDomains = 1;
    for (int i = 0; i < 8; ++i) {
        if (numThreads <= numThreadOptions[i] && numThreads > numThreadOptions[i - 1]) {
            numDomains = numThreadOptions[i];
        }
    }
    if (numThreads > 56) {
        numDomains = 56;
    }

    switch (numDomains) {
        case 1:
            return {1, 1, 1};
        case 2:
            return {2, 1, 1};
        case 4:
            return {2, 2, 1};
        case 8:
            return {4, 2, 1};
        case 14:
            return {7, 2, 1};
        case 16:
            return {4, 4, 1};
        case 28:
            return {7, 4, 1};
        case 56:
            return {8, 7, 1};
        default:
            return {1, 1, 1};

    }
}

void LinkedCellsContainer::initSubdomains(int numThreadsInit) {
    //now divide the domain into subdomains depending on the number of threads
    //  std::cout << "numThreads" << numThreadsInit << std::endl;
    std::array<double, 3> subdomainsPerDimension = computeSubdomainsPerDimension(numThreadsInit);
    Logger::logger->info("Subdomains per dimension: {} * {} * {}", subdomainsPerDimension[0], subdomainsPerDimension[1], subdomainsPerDimension[2]);
    //TODO: what if the amount of subdomains isn't a multiple of the amount of cells --> smaller subdomains
    //TODO: include the halo cells in the subdomain? --> yes

    // the cellsPerSubdomain variables denote the maximal amount of cells per subdomain
    // subdomains can have less cells if the following division don't yield an integral value
    auto cellsPerSubdomainX = std::ceil((domain_num_cells[0] + 2.0) / subdomainsPerDimension[0]);
    auto cellsPerSubdomainY = std::ceil((domain_num_cells[1] + 2.0) / subdomainsPerDimension[1]);
    auto cellsPerSubdomainZ = std::ceil((domain_num_cells[2] + 2.0) / subdomainsPerDimension[2]);
    auto numSubdomains = subdomainsPerDimension[0] * subdomainsPerDimension[1]
                         * subdomainsPerDimension[2];
    Logger::logger->info("cellsPerSubdomainX: {}", cellsPerSubdomainX);
    Logger::logger->info("cellsPerSubdomainY: {}", cellsPerSubdomainY);
    Logger::logger->info("cellsPerSubdomainZ: {}", cellsPerSubdomainZ);


    // instantiate numSubdomains subdomains
    for (unsigned int i = 0; i < subdomainsPerDimension[0]; ++i) {
        for (unsigned int j = 0; j < subdomainsPerDimension[1]; ++j) {
            for (unsigned int k = 0; k < subdomainsPerDimension[2]; ++k) {
                // std::array<unsigned int, 3> index{i, j, k};
                /*  subdomains.emplace(index, new Subdomain(delta_t, gravityConstant, cutoff_radius,
                                                          std::make_unique<LinkedCellsContainer>(
                                                                  *this))); //std::unique_ptr<LinkedCellsContainer>(this)));
                  subdomainsVector.emplace_back(subdomains.at(index));*/
                subdomainsVector.emplace_back(new Subdomain(delta_t, gravityConstant, cutoff_radius));
                // std::cout << "after subdomains vector init" << std::endl;
                auto index = subdomainsVector.size() - 1;
                auto lMax = std::min(cellsPerSubdomainX, ((domain_num_cells[0] + 2) - i * cellsPerSubdomainX));
                auto mMax = std::min(cellsPerSubdomainY, ((domain_num_cells[1] + 2) - j * cellsPerSubdomainY));
                auto nMax = std::min(cellsPerSubdomainZ, ((domain_num_cells[2] + 2) - k * cellsPerSubdomainZ));
                /*   std::cout << "lMax" << lMax << std::endl;
                   std::cout << "mMax" << mMax << std::endl;
                   std::cout << "nMax" << nMax << std::endl;*/
                for (int l = 0; l < lMax; ++l) {
                    for (int m = 0; m < mMax; ++m) {
                        for (int n = 0; n < nMax; ++n) {
                            // case distinction whether it is a cell at subdomain border or not
                            // for the addCell method
                            /*    std::cout << "cellCoordX: " << i * subdomainsPerDimension[0] + l << std::endl;
                                std::cout << "cellCoordY: " << j * subdomainsPerDimension[1] + m << std::endl;
                                std::cout << "cellCoordZ: " << k * subdomainsPerDimension[2] + n << std::endl;*/

                            //  std::cout<<"Index: "<<index<<"\n";
                            if (l == 0 || l == lMax - 1 || m == 0 || m == mMax - 1 ||
                                n == 0 || n == nMax - 1) {
                                int cellIndex = cellCoordToCellIndex(i * cellsPerSubdomainX + l - 1,
                                                                     j * cellsPerSubdomainY + m - 1,
                                                                     k * cellsPerSubdomainZ + n - 1);
                                //  std::cout<<"CellIndex: "<<cellIndex<<"\n";
                                subdomainsVector.at(index)->addCell(true, &cells.at(
                                        cellIndex));
                            } else {
                                subdomainsVector.at(index)->addCell(false, &cells.at(
                                        cellCoordToCellIndex(i * cellsPerSubdomainX + l - 1,
                                                             j * cellsPerSubdomainY + m - 1,
                                                             k * cellsPerSubdomainZ + n - 1)));
                            }
                        }
                    }
                }
            }
        }
    }
   /*for (int i = 0; i < subdomainsVector.size(); ++i) {
        std::cout << "Subdomain " << i << "consists of the following cells" << std::endl;
        for (auto &cell: subdomainsVector.at(i)->subdomainCells) {
            std::cout << "Cell ID: " << cell.second->getCellIndex() << std::endl;
        }
    }*/
   if(checkNoDuplicateCellsInDiffSubdomains()) {
       Logger::logger->info("There are no duplicate cells in the subdomains.");
   }
   else {
        Logger::logger->info("There are duplicate cells in the subdomains!");
   }
}

bool printCellNeighboursToComputeForcesWith() {

}

bool LinkedCellsContainer::checkNoDuplicateCellsInDiffSubdomains() {
    for (int i = 0; i < subdomainsVector.size() - 1; i++) {
        for (auto &cell: subdomainsVector.at(i)->subdomainCells) {
            int id_one = cell.second->getCellIndex();
            for (int j = i + 1; j < subdomainsVector.size(); ++j) {
                for (auto &cell2: subdomainsVector.at(j)->subdomainCells) {
                    int id_two = cell2.second->getCellIndex();
                    if (id_two == id_one) {
                        return false;
                    }
                }
            }
        }
    }
    return true;
}


void LinkedCellsContainer::initCellNeighbourReferences() {
    /*
     * addition: add the computeForcesWithIndices indices
     */
    // Loop through each cell according to their cell coordinates
    for (int cx = -1; cx < domain_num_cells[0] + 1; ++cx) {
        for (int cy = -1; cy < domain_num_cells[1] + 1; ++cy) {
            for (int cz = -1; cz < domain_num_cells[2] + 1; ++cz) {
                Cell &cell = cells.at(cellCoordToCellIndex(cx, cy, cz));

                // Loop through each of the current cells neighbour cells according to their cell coordinates
                // except the current cell itself
                for (int dx = -1; dx < 2; ++dx) {
                    for (int dy = -1; dy < 2; ++dy) {
                        for (int dz = -1; dz < 2; ++dz) {
                            if (dx == 0 && dy == 0 && dz == 0) continue;

                            // Get the cell index of the current neighbour cell
                            int cell_index = cellCoordToCellIndex(cx + dx, cy + dy, cz + dz);

                            // If the neighbour cell would be out of bounds, skip it
                            if (cell_index == -1) continue;

                            // Add the neighbour to the current cells neighbour references
                            // add to the current cell the indices of the neighbours which are
                            /*
                             *  add    add     add
                             *         cell    add
                             *
                             * + behind the cell (dz == 1)
                             * (needed for the subdomain parallization strategy)
                             */




                            if (dz == 1 || (dx == 1 && dy == 0 && dz == 0) || (dx == 1 && dy == 1 && dz == 0) || (dx == 0 && dy == 1 && dz == 0)
                                || (dx == -1 && dy == 1 && dz == 0)) {
                                cell.addToNeighboursToComputeForcesWith(&cells.at(cell_index));
                            }
                            Cell &curr_neighbour = cells.at(cell_index);
                            cell.addNeighbourReference(&curr_neighbour);
                        }
                    }
                }
            }
        }
    }
}

void LinkedCellsContainer::updateCellsParticleReferences() {
//#pragma omp parallel for schedule(dynamic)
    for (Cell &cell: cells) {
        cell.clearParticleReferences();
    }

    // clear the set of used cells
    occupied_cells_references.clear();
    occupied_cells_map.clear();

//#pragma omp parallel for schedule(dynamic)
    // add the particle references to the cells
    for (Particle &p: particles) {
        Cell *cell = particlePosToCell(p.getX());
        if (!occupied_cells_map.contains(cell->getCellIndex())) {
            occupied_cells_map.emplace(cell->getCellIndex(), cell);
        }
        occupied_cells_references.insert(cell);
        cell->addParticleReference(&p);
    }
}

void LinkedCellsContainer::deleteHaloParticles() {
    for (Cell *cell: halo_cell_references) {
        for (Particle *p: cell->getParticleReferences()) {
            particles.erase(std::find(particles.begin(), particles.end(), *p));
        }
    }
}

std::map<std::array<unsigned int, 3>, Subdomain *> LinkedCellsContainer::getSubdomains() {
    return subdomains;
}

std::vector<Subdomain *> LinkedCellsContainer::getSubdomainsVector() {
    return subdomainsVector;
}

void LinkedCellsContainer::parallel_step(const std::vector<std::shared_ptr<SimpleForceSource>> &simple_force_sources,
                                         const std::vector<std::shared_ptr<PairwiseForceSource>> &pairwise_force_sources,
                                         double delta_t, double gravityConstant, int strategy) {
#ifdef _OPENMP


    if (strategy == 1) { // parallelization strategy 1: subdomains
#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < getSubdomainsVector().size(); ++i) {
            for (auto &cell: getSubdomainsVector().at(i)->subdomainCells) {
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

        prepareForceCalculation();
#pragma omp for schedule(dynamic)
        for (auto *subdomain: subdomainsVector) {
            for (auto &cell: subdomain->subdomainCells) {
                if (cell.second->getCellType() != Cell::CellType::HALO) { // for particles in the halo cells no gravitational
                    //force has to be applied
                    for (auto *particle: cell.second->getParticleReferences()) { //directly apply gravitational
                        //force here for performance
                        auto currentF = particle->getF();
                        currentF[1] += particle->getM() * (-gravityConstant);
                        particle->setF(currentF);
                    }
                }
            }
        }

        ReflectiveBoundaryType::applyBoundaryConditions(*this);
        OutflowBoundaryType::applyBoundaryConditions(*this);
        PeriodicBoundaryType::applyBoundaryConditions(*this);

#pragma omp parallel for schedule(dynamic) //num_threads(8)
        for (auto *subdomain: subdomainsVector) {
            // std::cout << "num threads parallel_step: " << omp_get_num_threads() << std::endl;
            /*for (auto &cell: subdomain->subdomainCells) {
                if (cell.first) {
                    omp_set_lock(cell.second->getLock());
                }
                for (auto *particle: cell.second->getParticleReferences()) { //directly apply gravitational
                    //force here for performance
                    auto currentF = particle->getF();
                    currentF[1] += particle->getM() * (-gravityConstant);
                    particle->setF(currentF);
                }
                if (cell.first) {
                    omp_unset_lock(cell.second->getLock());
                }
            }*/
            // apply pairwise forces
            //TODO: maybe domain occupied_cell_references
            for (auto &cell: subdomain->subdomainCells) {
                /*
                only lock the cells which are at domain borders, because these are the only cells
                where race conflicts can occur(because every domain has at most one thread, so only one thread will
                work on the inner subdomain cells)
                whether the cell is at the subdomain border can be deduced from cell.first*/
                //

                //if (cell.first) {
                omp_set_lock(cell.second->getLock());
                //}
                if (cell.second->getCellType() != Cell::CellType::HALO) { // forces within the halo cell don't have to
                    //be calculated
                    for (auto it1 = cell.second->getParticleReferences().begin();
                         it1 != cell.second->getParticleReferences().end(); ++it1) {
                        Particle *p = *it1;
                        // calculate the forces between the particle and the particles in the same cell
                        // uses direct sum with newtons third law
                        for (auto it2 = (it1 + 1); it2 != cell.second->getParticleReferences().end(); ++it2) {
                            Particle *q = *it2;
                            std::array<double, 3> total_force{0, 0, 0};
                            for (auto &force: pairwise_force_sources) {  //maybe inline function for the force ?
                                total_force = total_force + force->calculateForce(*p, *q);
                            }
                            p->setF(p->getF() + total_force);
                            q->setF(q->getF() - total_force);
                        }
                    }
                }
                //if (cell.first) {
                omp_unset_lock(cell.second->getLock());
                //}
                // calculate the forces between the particles in the current cell
                // and particles in the neighbouring cells
                //  std::cout << "before neighbour force calc" << std::endl;
                for (Cell *neighbour: cell.second->getNeighboursToComputeForcesWith()) {
                    //if (neighbour.first) {
                    /// think a global lock order is established because of the deterministic force calculation
                    if (cell.second->getCellIndex() < neighbour->getCellIndex()) {
                        omp_set_lock(cell.second->getLock());
                        omp_set_lock(neighbour->getLock());
                    } else if (cell.second->getCellIndex() == neighbour->getCellIndex()) {
                        throw std::out_of_range("Different subdomains contain the same cell!");
                    } else {
                        omp_set_lock(neighbour->getLock());
                        omp_set_lock(cell.second->getLock());
                    }
                    //}
                    for (Particle *p: cell.second->getParticleReferences()) {
                        for (Particle *neighbour_particle: neighbour->getParticleReferences()) {
                            if (ArrayUtils::L2Norm(p->getX() - neighbour_particle->getX()) > cutoff_radius) continue;
                            for (const auto &force_source: pairwise_force_sources) {
                                std::array<double, 3> force = force_source->calculateForce(*p, *neighbour_particle);
                                p->setF(p->getF() + force);
                                neighbour_particle->setF(neighbour_particle->getF() - force);
                            }
                        }
                    }
                    if (cell.second->getCellIndex() < neighbour->getCellIndex()) {
                        omp_unset_lock(cell.second->getLock());
                        omp_unset_lock(neighbour->getLock());
                    } else {
                        omp_unset_lock(neighbour->getLock());
                        omp_unset_lock(cell.second->getLock());
                    }
                    // free the lock for the neighbour
                    // if (cell.first) {
                    /*   omp_unset_lock(cell.second->getLock());
                       omp_unset_lock(neighbour->getLock());*/
                }
            }
        }

        deleteHaloParticles();
        updateCellsParticleReferences();

#pragma omp parallel for schedule(dynamic)
        for (auto *subdomain: subdomainsVector) {
            for (auto &cell: subdomain->subdomainCells) {
                for (auto *particle: cell.second->getParticleReferences()) {
                    const std::array<double, 3> new_v = particle->getV() + (delta_t / (2 * particle->getM())) *
                                                                           (particle->getF() + particle->getOldF());
                    particle->setV(new_v);
                }
            }
        }

    } else if (strategy == 2) { // parallelization strategy 2 : particles parallelization
        /**
         * @brief: for parallelizaton strategy 2, each cell will be processed sequentielly,
         * within the cells threads are spawned to update the particles
         */
        /// update the positions
#pragma omp parallel for schedule(dynamic)
        for (auto &p: *this) {
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
            this->prepareForceCalculation();
        }

        /// apply simple forces
#pragma omp parallel for schedule(dynamic)
        for (auto &p: *this) { //directly apply gravitational
            //force here for performance
            auto currentF = p.getF();
            currentF[1] += p.getM() * (-gravityConstant);
            p.setF(currentF);
        }

        /// apply pairwise forces
        applyPairwiseForcesOptimized(pairwise_force_sources);
        //update the velocity
#pragma omp parallel for schedule(dynamic)
        for (auto &p: *this) {
            const std::array<double, 3> new_v = p.getV() + (delta_t / (2 * p.getM())) * (p.getF() + p.getOldF());
            p.setV(new_v);
        }
    } else if (strategy == 3) { //cell based parallelization, basic idea is that the occupied cells are
        // distributed on the threads

#pragma omp parallel for schedule(dynamic)
        for (auto &cell: cells) { // doesn't work as occupied_cell_references is an unordered_set
            for (auto *p: cell.getParticleReferences()) {
                const std::array<double, 3> new_x =
                        p->getX() + delta_t * p->getV() + (delta_t * delta_t / (2 * p->getM())) * p->getF();
                p->setX(new_x);

                // reset forces
                p->setOldF(p->getF());
                p->setF({0, 0, 0});
            }
        }
        prepareForceCalculation();

        ReflectiveBoundaryType::applyBoundaryConditions(*this);
        OutflowBoundaryType::applyBoundaryConditions(*this);
        PeriodicBoundaryType::applyBoundaryConditions(*this);

#pragma omp parallel for schedule(dynamic)
        for (auto &cell: cells) {
            omp_set_lock(cell.getLock());
            for (auto *particle: cell.getParticleReferences()) { //directly apply gravitational
                //force here for performance
                auto currentF = particle->getF();
                currentF[1] += particle->getM() * (-gravityConstant);
                particle->setF(currentF);
            }
            omp_unset_lock(cell.getLock());
        }

#pragma omp parallel for schedule(dynamic)
        for (auto &cell: cells) {

            /*
            only lock the cells which are at domain borders, because these are the only cells
            where race conflicts can occur(because every domain has at most one thread, so only one thread will
            work on the inner subdomain cells)
            whether the cell is at the subdomain border can be deduced from cell.first*/


            omp_set_lock(cell.getLock());
            for (auto it1 = cell.getParticleReferences().begin();
                 it1 != cell.getParticleReferences().end(); ++it1) {
                Particle *p = *it1;
                // calculate the forces between the particle and the particles in the same cell
                // uses direct sum with newtons third law
                //simple force

                /* auto currentF = p->getF();
                 currentF[1] += p->getM() * (-gravityConstant);
                 p->setF(currentF);*/

                for (auto it2 = (it1 + 1); it2 != cell.getParticleReferences().end(); ++it2) {
                    Particle *q = *it2;
                    std::array<double, 3> total_force{0, 0, 0};
                    for (auto &force: pairwise_force_sources) {  //maybe inline function for the force ?
                        total_force = total_force + force->calculateForce(*p, *q);
                    }
                    p->setF(p->getF() + total_force);
                    q->setF(q->getF() - total_force);
                }
            }
            omp_unset_lock(cell.getLock());
            // calculate the forces between the particles in the current cell
            // and particles in the neighbouring cells
            //  std::cout << "before neighbour force calc" << std::endl;
            for (Cell *neighbour: cell.getNeighboursToComputeForcesWith()) {
                //if (neighbour.first) {
                /// think a global lock order is established because of the deterministic force calculation
                if (cell.getCellIndex() < neighbour->getCellIndex()) {
                    omp_set_lock(cell.getLock());
                    omp_set_lock(neighbour->getLock());
                } else if (cell.getCellIndex() == neighbour->getCellIndex()) {
                    throw std::out_of_range("Different subdomains contain the same cell!");
                } else {
                    omp_set_lock(neighbour->getLock());
                    omp_set_lock(cell.getLock());
                }
                //}
                for (Particle *p: cell.getParticleReferences()) {
                    for (Particle *neighbour_particle: neighbour->getParticleReferences()) {
                        if (ArrayUtils::L2Norm(p->getX() - neighbour_particle->getX()) > cutoff_radius) continue;
                        for (const auto &force_source: pairwise_force_sources) {
                            std::array<double, 3> force = force_source->calculateForce(*p, *neighbour_particle);
                            p->setF(p->getF() + force);
                            neighbour_particle->setF(neighbour_particle->getF() - force);
                        }
                    }
                }
                if (cell.getCellIndex() < neighbour->getCellIndex()) {
                    omp_unset_lock(cell.getLock());
                    omp_unset_lock(neighbour->getLock());
                } else {
                    omp_unset_lock(neighbour->getLock());
                    omp_unset_lock(cell.getLock());
                }
            }
        }
        deleteHaloParticles();
        updateCellsParticleReferences();

#pragma omp for schedule(dynamic)
        for (auto &cell: cells) {
            for (auto &p: cell.getParticleReferences()) {
                const std::array<double, 3> new_v =
                        p->getV() + (delta_t / (2 * p->getM())) * (p->getF() + p->getOldF());
                p->setV(new_v);
            }
        }
    } else if (strategy == 4) {
#pragma omp parallel for schedule(dynamic)
        for (auto &cell: cells) { // doesn't work as occupied_cell_references is an unordered_set
            for (auto *p: cell.getParticleReferences()) {
                const std::array<double, 3> new_x =
                        p->getX() + delta_t * p->getV() + (delta_t * delta_t / (2 * p->getM())) * p->getF();
                p->setX(new_x);

                // reset forces
                p->setOldF(p->getF());
                p->setF({0, 0, 0});
            }
        }
        prepareForceCalculation();

        ReflectiveBoundaryType::applyBoundaryConditions(*this);
        OutflowBoundaryType::applyBoundaryConditions(*this);
        PeriodicBoundaryType::applyBoundaryConditions(*this);

#pragma omp parallel
        {
#pragma omp single
            {
                for (auto *cell: occupied_cells_references) {
#pragma omp task
                    processCell(cell);
                }
            }
        }

        deleteHaloParticles();
        updateCellsParticleReferences();

#pragma omp for schedule(dynamic)
        for (auto &cell: cells) {
            for (auto &p: cell.getParticleReferences()) {
                const std::array<double, 3> new_v =
                        p->getV() + (delta_t / (2 * p->getM())) * (p->getF() + p->getOldF());
                p->setV(new_v);
            }
        }

    }
#endif
}

void LinkedCellsContainer::processCell(Cell *cell) {
#ifdef _OPENMP
    omp_set_lock(cell->getLock());
    for (auto it1 = cell->getParticleReferences().begin();
         it1 != cell->getParticleReferences().end(); ++it1) {
        Particle *p = *it1;
        // calculate the forces between the particle and the particles in the same cell
        // uses direct sum with newtons third law
        //simple force
        auto currentF = p->getF();
        currentF[1] += p->getM() * (-gravityConstant);
        p->setF(currentF);


        for (auto it2 = (it1 + 1); it2 != cell->getParticleReferences().end(); ++it2) {
            Particle *q = *it2;
            std::array<double, 3> total_force{0, 0, 0};
            for (auto &force: pairwise_force_sources) {  //maybe inline function for the force ?
                total_force = total_force + force->calculateForce(*p, *q);
            }
            p->setF(p->getF() + total_force);
            q->setF(q->getF() - total_force);
        }
    }
    omp_unset_lock(cell->getLock());
    // calculate the forces between the particles in the current cell
    // and particles in the neighbouring cells
    //  std::cout << "before neighbour force calc" << std::endl;
    for (Cell *neighbour: cell->getNeighboursToComputeForcesWith()) {
        //if (neighbour.first) {
        /// think a global lock order is established because of the deterministic force calculation
        if (cell->getCellIndex() < neighbour->getCellIndex()) {
            omp_set_lock(cell->getLock());
            omp_set_lock(neighbour->getLock());
        } else if (cell->getCellIndex() == neighbour->getCellIndex()) {
            throw std::out_of_range("Different subdomains contain the same cell!");
        } else {
            omp_set_lock(neighbour->getLock());
            omp_set_lock(cell->getLock());
        }
        for (Particle *p: cell->getParticleReferences()) {
            for (Particle *neighbour_particle: neighbour->getParticleReferences()) {
                if (ArrayUtils::L2Norm(p->getX() - neighbour_particle->getX()) > cutoff_radius) continue;
                for (const auto &force_source: pairwise_force_sources) {
                    std::array<double, 3> force = force_source->calculateForce(*p, *neighbour_particle);
                    p->setF(p->getF() + force);
                    neighbour_particle->setF(neighbour_particle->getF() - force);
                }
            }
        }
        if (cell->getCellIndex() < neighbour->getCellIndex()) {
            omp_unset_lock(cell->getLock());
            omp_unset_lock(neighbour->getLock());
        } else {
            omp_unset_lock(neighbour->getLock());
            omp_unset_lock(cell->getLock());
        }
    }
#endif
}

void LinkedCellsContainer::setPairwise(std::vector<std::shared_ptr<PairwiseForceSource>> pairwiseForceSources) {

}
