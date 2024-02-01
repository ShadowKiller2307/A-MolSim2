#pragma once

#include "physics/pairwiseforces/PairwiseForceSource.h"

/**
 *  @brief Class to calculate the harmonic force between particles. Implements the interface PairwiseForceSource.
 *
 * Implementation of the force calculation to simulate harmonic forces between particles.
 */
class HarmonicPotential : public SimpleForceSource
{
private:
	constexpr static const double sqrt2 = 1.4142135623730951;
	double r0;	   // equilibrium distance for orthaogonal particles
	double r0sqrt; // equilibrium distance for diagonal particles
	double k;	   // stiffness

public:
	/**
	 * @brief Constructor for the HarmonicPotential class
	 * @param r0 The equilibrium distance for orthaogonal particles
	 * @param k The stiffness
	 */
	explicit HarmonicPotential(double r0, double k);

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
	std::array<double, 3UL> calculateForce(Particle &p) const override;

	/**
	 * @brief Returns "Harmonic" as name of the force
	 */
	explicit
	operator std::string() const override;

	/**
	 * @brief Sets the equilibrium distance and stiffness
	 *
	 * @param r0 The equilibrium distance for orthaogonal particles
	 * @param k The stiffness
	 */
	void setR0andK(double r0, double k);
};