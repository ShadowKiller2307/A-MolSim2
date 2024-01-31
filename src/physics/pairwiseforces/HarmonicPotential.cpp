#include "HarmonicPotential.h"

#include "utils/ArrayUtils.h"

HarmonicPotential::HarmonicPotential(double r0, double k) : r0(r0), r0sqrt(r0 * sqrt2), k(k)
{
}

std::array<double, 3UL> HarmonicPotential::calculateForce(Particle &p, Particle &q) const { return {0.0, 0.0, 0.0}; };

std::array<double, 3UL> HarmonicPotential::calculateForce(Particle &p, Particle &q, bool is_diagonal) const
{
	const double eqilDist = is_diagonal ? r0sqrt : r0;
	const double distance = ArrayUtils::L2Norm(p.getX() - q.getX());

	return k * (distance - r0) * (1 / distance) * (p.getX() - q.getX());
}

HarmonicPotential::operator std::string() const { return "Harmonic"; }

void HarmonicPotential::setR0andK(double r0, double k)
{
	this->r0 = r0;
	this->r0sqrt = r0 * sqrt2;
	this->k = k;
};