#include "Config.hpp"

void Config::addStock(const Stock &stock) { stocks.push_back(stock); }

void Config::addProcess(const Process &process) {
    processes.push_back(process);
}

void Config::setOptimizeGoal(const std::vector<std::string> &goal) {
    optimizeGoal = goal;
}

const std::vector<Stock> &Config::getStocks() const { return stocks; }

const std::vector<Process> &Config::getProcesses() const { return processes; }

const std::vector<std::string> &Config::getOptimizeGoal() const {
    return optimizeGoal;
}