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

    /**
     * @brief: the pair consists of:
     * bool: is the cell located at the subdomain border
     * Cell*: pointer to the cell belonging to the subdomain
     */
    std::set<std::pair<bool, Cell*>> subdomainCells;

    /**
     * @brief initialize a subdomain
     * based on the number of threads subdomains are intialized which consist of multiple cells
     */
    Subdomain(double delta_t, double gravityConstant, double cutoffRadius);

    /**
     * @brief add a cell to the subdomain
     */
   void addCell(bool isAtSubdomainBorder, Cell* cellToAdd);

  /**
   * @brief updates the particle positions within the current domain
   * not used anymore as there were problems calling it in a pragma omp parallel section(now subdomains are
   * processed in parallel_step)
   */
  void updateParticlePositions();
  /**
   * @updates the subdomain by applying the remaining operations needed on the particles
   * -> apply simple forces
   * -> apply pairwise forces
   * -> update the velocity
   * also not used anymore
   */
  void updateSubdomain(const std::vector<std::shared_ptr<PairwiseForceSource>>& force_sources);

  void calculateForcesBetweenCells(Cell* one, Cell* two);

private:
  //  std::set<std::pair<bool, Cell*>> subdomainCells; // a set of the subdomains
    double delta_t;
    double gravityConstant;
    double cutoffRadius;
};

