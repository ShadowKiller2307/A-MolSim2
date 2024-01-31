#pragma once
#include <omp.h>
#include <set>
#include "particles/containers/linkedcells/cells/Cell.h"
#include <omp.h>
#include <memory>
#include "utils/ArrayUtils.h"
#include "physics/pairwiseforces/PairwiseForceSource.h"

// class for the subdomains, which are part of parallelization strategy 1
//template<unsigned N>
//class LinkedCellsContainer;

class Subdomain{
public:

    std::set<std::pair<bool, Cell*>> subdomainCells;

    Subdomain(double delta_t, double gravityConstant, double cutoffRadius);
    /**
     * @brief: based on the number of threads subdomains are intialized which consist of multiple cells
     */
  //  void initializeSubdomains();
  void addCell(bool isAtSubdomainBorder, Cell* cellToAdd);

  /**
   * @brief updates the particle positions within the current domain
   */
  void updateParticlePositions();
  /**
   * @updates the subdomain by applying the remaining operations needed on the particles
   * -> apply simple forces
   * -> apply pairwise forces
   * -> update the velocity
   */
  void updateSubdomain(const std::vector<std::shared_ptr<PairwiseForceSource>>& force_sources);

  void calculateForcesBetweenCells(Cell* one, Cell* two);

private:
    /**
     * @brief: the pair consists of:
     * bool: is the cell located at the subdomain border
     * Cell*: pointer to the cell belonging to the subdomain
     */
  //  std::set<std::pair<bool, Cell*>> subdomainCells; // a set of the subdomains
    double delta_t;
    double gravityConstant;
    double cutoffRadius;
   /* std::set<std::pair<bool, Cell*>> subdomainCells;*/
  /*  std::array<double, 3> domainSize;
    std::array<unsigned, 3> cellsPerDimension;
    double cutoffRadius;
    unsigned numThreads;*/

  // friend LinkedCellsContainer;
};

