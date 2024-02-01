#include "PullUpForce.h"

#include "utils/ArrayUtils.h"

PullUpForce::PullUpForce(double intensity) : intensity(intensity) {}

std::array<double, 3UL> PullUpForce::calculateForce(Particle &p) const
{
	if (p.getType() == 1)
	{
		return {0, 0, p.getM() * intensity};
	}
	return {0, 0, 0};
}

PullUpForce::operator std::string() const { return "PullUp"; }
