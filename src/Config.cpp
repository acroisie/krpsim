#include "Config.hpp"

void Config::addStock(const Stock &stock) {
    stocks_.push_back(stock);
}

void Config::addProcess(const Process &process) {
    processes_.push_back(process);
}

void Config::setOptimizeGoal(const std::vector<std::string> &goal) {
    optimizeGoal_ = goal;
}

const std::vector<Stock> &Config::getStocks() const {
    return stocks_;
}

const std::vector<Process> &Config::getProcesses() const {
    return processes_;
}

const std::vector<std::string> &Config::getOptimizeGoal() const {
    return optimizeGoal_;
}
