#include "HarmonicPotential.h"

#include "utils/ArrayUtils.h"

std::array<double, 3UL> HarmonicPotential::calculateForce(Particle &p, Particle &q) const { return {0.0, 0.0, 0.0}; };

std::array<double, 3UL> HarmonicPotential::calculateForce(Particle &p, Particle &q, bool is_diagonal) const
{
	// TODO: Hardcoded
	const double r0 = is_diagonal ? 2.2 * sqrt2 : 2.2; // equilibrium distance
	const double k = 300;							   // stiffness
	const double distance = ArrayUtils::L2Norm(p.getX() - q.getX());
	return k * (distance - r0) * (1 / distance) * (p.getX() - q.getX());
}

HarmonicPotential::operator std::string() const { return "Harmonic"; };