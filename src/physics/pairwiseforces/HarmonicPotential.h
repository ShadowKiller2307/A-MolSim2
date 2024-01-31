#pragma once

#include "physics/pairwiseforces/PairwiseForceSource.h"

/**
 *  @brief Class to calculate the harmonic force between particles. Implements the interface PairwiseForceSource.
 *
 * Implementation of the force calculation to simulate harmonic forces between particles.
 */
class HarmonicPotential : public PairwiseForceSource
{
private:
	constexpr static const double sqrt2 = 1.4142135623730951;

public:
	std::array<double, 3UL> calculateForce(Particle &p, Particle &q) const override;

	/**
	 * @brief Calculates the harmonic forces between two orthogonal or diagonal particles
	 *
	 * @param p Particle
	 * @param q Particle
	 * @param is_diagonal true if p and q are diagonal, false if they are orthogonal
	 * @return harmonic force exerted by q on p
	 *
	 * Calculates the harmonic force which q exerts on p
	 */
	std::array<double, 3UL> calculateForce(Particle &p, Particle &q, bool is_diagonal) const;

	/**
	 * @brief Returns "Harmonic" as name of the force
	 */
	explicit
	operator std::string() const override;
};