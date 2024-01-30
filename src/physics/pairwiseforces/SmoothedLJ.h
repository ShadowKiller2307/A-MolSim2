#pragma once

#include "physics/pairwiseforces/PairwiseForceSource.h"

class SmoothedLJ : public PairwiseForceSource{
public:
    SmoothedLJ(double rl, double cutOffRadius);

    /**
     *@brief calculates the smoothed LJ force
     *
     * @param p first particle for the force calculation
     * @param q second particle for the force calculation
     * @return The smoothed Lennard Jones Forces that acts between the two particles
     */
    std::array<double, 3UL> calculateForce(Particle &p, Particle& q) const override;

    /**
     * @brief Returning "Smoothed LJ" as the name of the force
     * @return "Smoothed LJ"
     */
    explicit operator std::string() const override;

private:
    double _r_l;
    double _cutOffRadius;
};

