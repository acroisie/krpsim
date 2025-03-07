#pragma once
#include <vector>
#include <string>
#include <map>
#include "Config.hpp"

class ProcessManager {
public:
    ProcessManager(const Config &config, int delayLimit);
    bool runSimulation();
    bool runGeneticAlgorithm();

private:
    const Config &config_;
    int delayLimit_;
    int currentCycle_;

    std::map<std::string, int> currentStocks_;
    std::vector<std::pair<int, std::string>> executionLogs_;

    std::vector<const Process*> getRunnableProcesses();
    bool executeProcess(const Process* process);
    void updateStocksWithOutputs();
    void generateOutput();
};
