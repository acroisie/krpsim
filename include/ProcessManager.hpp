#pragma once

#include "Config.hpp"
#include "Simulator.hpp"
#include "Optimizer.hpp"

class ProcessManager {
public:
    ProcessManager(const Config &config, int delayLimit);
    void runGeneticAlgorithm();

private:
    Config config_;
    Simulator simulator_;
    Optimizer optimizer_;
};