#include "DiffusionInterceptor.h"
#include "utils/ArrayUtils.h"
#include "simulation/SimulationParams.h"
DiffusionInterceptor::DiffusionInterceptor() {
    test_mode = false;
};

void DiffusionInterceptor::onSimulationStart(Simulation &simulation) {

    csv_writer = std::make_unique<CSVWriter>(simulation.params.output_dir_path / "diffusion.csv");
    csv_writer->initialize({"iteration","variance"});
    every_nth_iteration = 1000;



    for(auto& pos:simulation.initial_pos_of_particles){
        previous_references.push_back(pos);
    }

}

void DiffusionInterceptor::operator()(size_t iteration, Simulation& simulation) {

    double var = 0;
    auto number_of_particles = static_cast<double>(simulation.particle_container->size());
    double coefficient = 1.0/number_of_particles;

    double sum = 0;
    size_t i = 0;
    for (auto & particle : *simulation.particle_container) {
        std::array<double,3> displacement = particle.getX()-previous_references.at(i);
        std::cout<<"Previous reference: "<<previous_references.at(i)<<std::endl;
        std::cout<<"Displacement to add: "<<particle.getDisplacementToAdd()<<std::endl;

        if(test_mode){
            displacement.at(0)+=6;
        }
        else{
            displacement= displacement+particle.getDisplacementToAdd();
        }

        particle.setDisplacementToAdd({0,0,0});
        double norm = ArrayUtils::L2Norm(displacement);
        norm *=norm;


        previous_references.at(i) = particle.getX();
        sum+=norm;
        i++;
    }
    var = coefficient*sum;
    Logger::logger->error("coeff: {}, sum: {}\n",coefficient,sum);
    if(test_mode){
        variances_for_test_mode.push_back(var);
    }

    csv_writer->writeRow({iteration,var});

}

void DiffusionInterceptor::onSimulationEnd(size_t iteration, Simulation &simulation) {}

DiffusionInterceptor::operator std::string() const {
    return "Diffusion";
}

void DiffusionInterceptor::logSummary(int depth) const {
    std::string indent = std::string(depth * 2, ' ');

    Logger::logger->info("{}╟┤{}Diffusion: {}", indent, ansi_orange_bold, ansi_end);
}

void DiffusionInterceptor::setTestMode(bool mode) {
    test_mode = mode;
    if(mode) {
        variances_for_test_mode = std::vector<double>();
        every_nth_iteration = 7;
    }
}

std::vector<double> DiffusionInterceptor::getVariancesForTestMode() const {
    return variances_for_test_mode;
}

bool DiffusionInterceptor::isTestMode() const {
    return test_mode;
}

void DiffusionInterceptor::setPreviousReferences(std::vector<std::array<double, 3>>& ref) {
    if(test_mode) {
        previous_references = ref;
    }
}


