#include "SmoothedLJ.h"
#include "utils/ArrayUtils.h"

std::array<double, 3UL> SmoothedLJ::calculateForce(Particle &p, Particle &q) const
{
    // based on the distance between particles p and j a different formula
    // for the forces will be applied
    double cutoffRadius = _cutOffRadius;
    double r_l = _r_l;

    const auto displacement = q.getX() - p.getX();
    const double distancePQ = ArrayUtils::L2Norm(displacement);

    const double sigma = (p.getSigma() + q.getSigma()) / 2;
    const double epsilon = std::sqrt(p.getEpsilon() * q.getEpsilon());
    // case d_ij <= r_l

    if (distancePQ <= r_l)
    {
        // old lennard jones force calculation
        const auto f_lennard_jones =
            (24 * epsilon / (std::pow(distancePQ, 2))) * (std::pow(sigma / distancePQ, 6) - 2 * std::pow(sigma / distancePQ, 12)) * displacement;
        return f_lennard_jones;
    }
    // case r_l <= d_ij <= r_c
    else
    {
        const auto new_part =
            (-((24 * std::pow(sigma, 6) * epsilon) / (std::pow(distancePQ, 14) / std::pow(cutoffRadius - r_l, 3))) * (cutoffRadius - distancePQ)) * displacement;
        auto summandOne = std::pow(cutoffRadius, 2) * (2 * std::pow(sigma, 6) - std::pow(distancePQ, 6));
        auto summandTwo = cutoffRadius * (3 * r_l - distancePQ) * (std::pow(distancePQ, 6) - 2 * std::pow(sigma, 6));
        auto summandThree = distancePQ * (5 * r_l * std::pow(sigma, 6) - 2 * r_l * std::pow(distancePQ, 6) - 3 * std::pow(sigma, 6) * distancePQ + std::pow(distancePQ, 7));
        return (summandOne + summandTwo + summandThree) * new_part;
    }
    // case d_ij > r_c doesn't need to be considered, as the force will never
    // be calculated for particles with a distance greater than the cutoff radius
}

SmoothedLJ::operator std::string() const
{
    return "Smoothed LJ";
}

SmoothedLJ::SmoothedLJ(double rl,double cutOffRadius) {
    _r_l = rl;
    _cutOffRadius = cutOffRadius;
}
