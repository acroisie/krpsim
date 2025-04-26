#pragma once
#include "Config.hpp"
#include "GeneticAlgorithm.hpp"
#include "Individual.hpp"
#include "Simulator.hpp"
#include <map>
#include <string>
#include <vector>

class ProcessManager {
  public:
    ProcessManager(const Config &config, int timeLimit);
    void run();

  private:
    const Config &config;
    Simulator simulator;
    GeneticAlgorithm geneticAlgorithm;
    int timeLimit;

    static int calculateMaxSequenceLength(const Config &config, int timeLimit);

    void generateOutput(const Individual &bestSolution);
};