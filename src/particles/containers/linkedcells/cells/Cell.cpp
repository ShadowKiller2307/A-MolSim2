#include "Cell.h"

Cell::Cell(CellType cell_type) : cell_type(cell_type) {
    neighbour_references.reserve(26);
    already_influenced_by.reserve(16);
    omp_init_lock(&cellLock);
}

Cell::CellType Cell::getCellType() const { return cell_type; }

std::vector<Particle*>& Cell::getParticleReferences() { return particle_references; }

std::vector<Cell*>& Cell::getNeighbourReferences() { return neighbour_references; }

std::unordered_set<Cell*>& Cell::getAlreadyInfluencedBy() { return already_influenced_by; }

void Cell::addParticleReference(Particle* p) { particle_references.push_back(p); }

void Cell::clearParticleReferences() { particle_references.clear(); }

void Cell::addNeighbourReference(Cell* c) { neighbour_references.push_back(c); }

void Cell::addAlreadyInfluencedBy(Cell* c) { already_influenced_by.insert(c); }

void Cell::clearAlreadyInfluencedBy() { already_influenced_by.clear(); }

void Cell::addToNeighboursToComputeForcesWith(Cell* cell) {
    neighboursToComputeForcesWith.emplace(cell);
}

std::unordered_set<Cell*> Cell::getNeighboursToComputeForcesWith() {
    return neighboursToComputeForcesWith;
}

Cell::~Cell() {
    omp_destroy_lock(&cellLock);
}

omp_lock_t* Cell::getLock() {
    return &cellLock;
}

void Cell::calculateForcesBetweenCells(Cell *other) {


}
