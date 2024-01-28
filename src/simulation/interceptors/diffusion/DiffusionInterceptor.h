#pragma once
#include <chrono>
#include <memory>
#include <iostream>

#include "io/csv/CSVWriter.h"
#include "simulation/interceptors/SimulationInterceptor.h"

class DiffusionInterceptor: public SimulationInterceptor{
public:
    DiffusionInterceptor();

    /***
     * Initializes the writer and gets the initial position of the particles from the simulation parameter
     * @param simulation Uses it to obtain the initial positions of the particles
     */
    void onSimulationStart(Simulation& simulation) override;

    /***
     *
     * @param iteration Calculation at a specific time = Calculation at a specific iteration
     * @param simulation Used for getting the particles
     */
    void operator()(size_t iteration, Simulation& simulation) override;

    /***
     * @brief Nothing to do here. Therefore it is empty.
     */
    void onSimulationEnd(size_t iteration, Simulation& simulation) override;
    /***
     * @brief Returns the name of the interceptor
     * @return Name of the interceptor
     */
    explicit operator std::string() const override;
    /***
     * @brief Logs the summary of the interceptor
     * @param depth Depth of the summary
     */
    void logSummary(int depth) const override;

    void setTestMode(bool mode);

    void setPreviousReferences(std::vector<std::array<double,3>>& references);

    [[nodiscard]] std::vector<double> getVariancesForTestMode() const;
    [[nodiscard]] bool isTestMode() const;

private:
    std::unique_ptr<CSVWriter> csv_writer;

    ///@brief Each time the operator is called the particle positions used in the calculations become the previous references
    std::vector<std::array<double,3>> previous_references;

    bool test_mode;
    std::vector<double> variances_for_test_mode;



};