#pragma once
#include "Config.hpp"
#include "Individual.hpp"
#include "Simulator.hpp"
#include "GeneticAlgorithm.hpp"
#include <map>
#include <string>
#include <vector>

class ProcessManager {
  public:
    ProcessManager(const Config &config, int delayLimit);
    void run();

  private:
    const Config &config_;
    Simulator simulator_;
    GeneticAlgorithm geneticAlgorithm_;
    
    std::map<std::string, int> currentStocks_;
    std::vector<std::pair<int, std::string>> executionLogs_;
    
    void generateOutput(const Individual& bestSolution);
};