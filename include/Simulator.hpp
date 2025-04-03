#pragma once

#include "Config.hpp"
#include "Individual.hpp"
#include <limits>
#include <map>
#include <queue>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

class Simulator {
  public:
    struct RunningProcess {
        const Process *processPtr = nullptr;
        int completionTime = 0;

        bool operator<(const RunningProcess &other) const {
            return completionTime > other.completionTime;
        }
    };

    struct SimulationResult {
        double fitnessScore = std::numeric_limits<double>::lowest();
        std::vector<std::pair<int, string>> executionLog;
        std::map<std::string, int> finalStocks;
        int finalCycle = 0;
        bool timeoutReached = false;

        SimulationResult() = default;
    };

    explicit Simulator(const Config &config, int timeLimit);

    SimulationResult runSimulation(const std::vector<std::string> &sequence);

  private:
    const Config &config;
    int timeLimit;
    std::unordered_map<std::string, 
};