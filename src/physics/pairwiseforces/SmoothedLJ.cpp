#include "SmoothedLJ.h"
#include "utils/ArrayUtils.h"

std::array<double, 3UL> calculateForce(Particle &p, Particle& q) const {
    // based on the distance between particles p and j a different formula
    // for the forces will be applied
    double cutoffRadius = 2.3; //TODO change later
    double r_l = 1.9;

    const auto displacement = q.getX() - p.getX();
    const double distancePQ = ArrayUtils::L2Norm(displacement);
    // TODO: maybe apply the mixing rule only if sigma and epsilon are different?
    const double sigma = (p.getSigma() + q.getSigma()) / 2;
    const double epsilon = std::sqrt(p.getEpsilon() * q.getEpsilon());
    // case d_ij <= r_l
    // TODO: maybe some parts of the calculation can be optimized(for example that every same calculation is only executed once)
    if (distancePQ <= r_l) {
        // old lennard jones force calculation
        const auto f_lennard_jones =
                (24 * epsilon / (std::pow(distancePQ, 2))) * (std::pow(sigma / distancePQ, 6) - 2 * std::pow(sigma / distancePQ, 12)) * displacement;

        const auto new_part =
                displacement * -((24* std::pow(sigma, 6) * epsilon)/(std::pow(distancePQ, 14) / std::pow(cutoffRadius - r_l, 3))) * (cutoffRadius - distancePQ);

        return f_lennard_jones * new_part;
    }
    // case r_l <= d_ij <= r_c
    else {
        auto summandOne = std::pow(cutoffRadius, 2) * (2 * std::pow(sigma, 6) - std::pow(distancePQ, 6));
        auto summandTwo = cutoffRadius * (3* r_l - distancePQ) * (std::pow(distancePQ, 6) - 2 * std::pow(sigma, 6));
        auto summandThree = distancePQ * (5*r_l*std::pow(sigma, 6) - 2 * r_l * std::pow(distancePQ, 6) - 3 * std::pow(sigma, 6) * distancePQ + std::pow(distancePQ, 7));
        return summandOne+summandTwo+summandThree;
    }
    // case d_ij > r_c doesn't need to be considered, as the force will never
    // be calculated for particles with a distance greater than the cutoff radius
}

explicit operator std::string() const {
    return "Smoothed LJ";
}