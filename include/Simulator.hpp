#pragma once
#include "Config.hpp"
#include "Process.hpp"
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

class Simulator {
  public:
    struct RunningProcess {
        const Process *processPtr;
        int completionTime;
        bool operator>(const RunningProcess &other) const {
            return completionTime > other.completionTime;
        }
    };

    struct SimulationResult {
        double fitness = std::numeric_limits<double>::lowest();
        std::vector<std::pair<int, std::string>> executionLog;
        std::map<std::string, int> finalStocks;
        int finalCycle = 0;
        bool timeoutReached = false;
    };

    Simulator(const Config &config, int timeLimit);

    SimulationResult runSimulation(const std::vector<std::string> &processSequence);

    const std::unordered_map<std::string, int> &getProcessPriority() const {
        return processPriority;
    }

  private:
    const Config &config;
    int timeLimit;

    std::unordered_map<std::string, const Process *> processMap;
    std::unordered_map<std::string, int> processPriority;

    void buildProcessPriority();
    bool canStartProcess(const Process *process, const std::map<std::string, int> &stocks) const;
    const Process *getProcessByName(const std::string &name) const;
};
