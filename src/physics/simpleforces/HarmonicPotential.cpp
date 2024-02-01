#include "HarmonicPotential.h"

#include "utils/ArrayUtils.h"

HarmonicPotential::HarmonicPotential(double r0, double k) : r0(r0), r0sqrt(r0 * sqrt2), k(k) {}

std::array<double, 3UL> HarmonicPotential::calculateForce(Particle &p) const
{
	std::array<double, 3> force = {0.0, 0.0, 0.0};
	for (auto &neighbour : p.getNeighbours())
	{
		const std::ptrdiff_t diff = neighbour.first;
		const bool is_diagonal = neighbour.second;
		const Particle &q = *(&p + diff);
		const double distance = ArrayUtils::L2Norm(p.getX() - q.getX());
		force = force + k * (distance - (is_diagonal ? r0sqrt : r0)) * (1 / distance) * (p.getX() - q.getX());
	}

	return force;
}

HarmonicPotential::operator std::string() const { return "Harmonic"; }

void HarmonicPotential::setR0andK(double r0, double k)
{
	this->r0 = r0;
	this->r0sqrt = r0 * sqrt2;
	this->k = k;
};