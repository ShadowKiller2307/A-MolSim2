#include "diffusion.h"
#include "utils/ArrayUtils.h"
DiffusionInterceptor::DiffusionInterceptor() {}

void DiffusionInterceptor::onSimulationStart(Simulation &simulation) {

    csv_writer = std::make_unique<CSVWriter>(simulation.params.output_dir_path / "diffusion.csv");
    csv_writer->initialize({"iteration","variance"});

    for(auto& pos:simulation.initial_pos_of_particles){
        previous_references.push_back(pos);
    }

}

void DiffusionInterceptor::operator()(size_t iteration, Simulation& simulation) {

    double var = 0;
    double coefficient = 1/simulation.particle_container->size();

    double sum = 0;
    size_t i = 0;
    for (auto it = simulation.particle_container->begin(); it != simulation.particle_container->end(); it++) {
        auto& particle = *it;
        double norm = ArrayUtils::L2Norm(particle.getX()-previous_reference.at(i));
        norm = std::pow(norm,1/2);
        previous_references.at(i) = particle.getX();
        sum+=norm;
        i++;

    }
    var = coefficient*sum;

    csv_writer->writeRow({iteration,var});

}

void DiffusionInterceptor::onSimulationEnd(size_t iteration, Simulation &simulation) {}