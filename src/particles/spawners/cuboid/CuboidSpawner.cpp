#include "CuboidSpawner.h"

#include "particles/Particle.h"
#include "physics/thermostats/Thermostat.h"
#include "utils/ArrayUtils.h"

CuboidSpawner::CuboidSpawner(const std::array<double, 3> &lower_left_corner, const std::array<int, 3> &grid_dimensions, double grid_spacing,
                             double mass, const std::array<double, 3> &initial_velocity, int type, double epsilon, double sigma,
                             bool third_dimension, double initial_temperature)
    : lower_left_corner(lower_left_corner),
      grid_dimensions(grid_dimensions),
      grid_spacing(grid_spacing),
      mass(mass),
      type(type),
      epsilon(epsilon),
      sigma(sigma),
      initial_velocity(initial_velocity),
      third_dimension(third_dimension),
      initial_temperature(initial_temperature) {}

int CuboidSpawner::spawnParticles(std::vector<Particle> &particles) const
{
    int offset = particles.size();
    particles.reserve(particles.size() + getEstimatedNumberOfParticles());
    for (int i = 0; i < grid_dimensions[0]; i++)
    {
        for (int j = 0; j < grid_dimensions[1]; j++)
        {
            for (int k = 0; k < grid_dimensions[2]; k++)
            {
                const auto grid_pos = std::array<double, 3>{static_cast<double>(i), static_cast<double>(j), static_cast<double>(k)};

                const auto x = lower_left_corner + grid_spacing * grid_pos;
                int actual_type = type;
                if (i > 16 && i < 19 && j > 23 && j < 26)
                {
                    actual_type = 1;
                }
                Particle particle(x, initial_velocity, mass, actual_type, epsilon, sigma);
                Thermostat::setParticleTemperature(initial_temperature, particle, third_dimension ? 3 : 2);
                particles.push_back(std::move(particle));
            }
        }
    }

    // set up neighbours in the membrane
    for (int i = 0; i < grid_dimensions[0]; i++)
    {
        for (int j = 0; j < grid_dimensions[1]; j++)
        {
            for (int k = 0; k < grid_dimensions[2]; k++)
            {
                const int indexP = getParticleIndexByPosition(offset, {i, j, k});
                // connect top three particles
                for (int l = -1; l < 2; l++)
                {
                    const int indexQ = getParticleIndexByPosition(offset, {i + l, j + 1, k});
                    if (indexQ == -1)
                    {
                        continue;
                    }
                    std::ptrdiff_t diff = &(particles[indexQ]) - &(particles[indexP]);
                    particles[indexP].addNeighbour(diff, l == 0 ? false : true);
                    particles[indexQ].addNeighbour(-diff, l == 0 ? false : true);
                }
                // connect right particle
                const int indexQ = getParticleIndexByPosition(offset, {i + 1, j, k});
                if (indexQ == -1)
                {
                    continue;
                }
                std::ptrdiff_t diff = &particles[indexQ] - &particles[indexP];
                particles[indexP].addNeighbour(diff, false);
                particles[indexQ].addNeighbour(-diff, false);
            }
        }
    }
    return grid_dimensions[0] * grid_dimensions[1] * grid_dimensions[2];
}

const int CuboidSpawner::getParticleIndexByPosition(const int offset, const std::array<double, 3> &position) const
{
    if (position[0] < 0 || position[0] >= grid_dimensions[0] || position[1] < 0 || position[1] >= grid_dimensions[1] || position[2] < 0 || position[2] >= grid_dimensions[2])
    {
        return -1;
    }
    return offset + grid_dimensions[2] * grid_dimensions[1] * position[0] + grid_dimensions[2] * position[1] + position[2];
}

size_t CuboidSpawner::getEstimatedNumberOfParticles() const
{
    return static_cast<size_t>(grid_dimensions[0]) * grid_dimensions[1] * grid_dimensions[2];
}