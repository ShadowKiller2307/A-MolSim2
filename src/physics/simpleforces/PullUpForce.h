#pragma once

#include "physics/simpleforces/SimpleForceSource.h"

/**
 *  @brief Class to pull up particles. Implements the interface SimpleForceSource.
 *
 * Implementation of the force calculation to simulate harmonic forces between particles.
 */
class PullUpForce : public SimpleForceSource
{
private:
	double intensity;

public:
	/**
	 * @brief Constructor for the PullUpForce class
	 */
	explicit PullUpForce(double intensity);

	/**
	 * @brief Pulls up the particle with type 1 by a given intensity
	 *
	 * @param p Particle
	 * @return force p is pulled up by
	 *
	 */
	std::array<double, 3UL> calculateForce(Particle &p) const override;

	/**
	 * @brief Returns "PullUp" as name of the force
	 */
	explicit
	operator std::string() const override;
};