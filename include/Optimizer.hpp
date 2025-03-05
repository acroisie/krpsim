#pragma once
#include <string>
#include <vector>
#include "Config.hpp"

class Optimizer {
	explicit Optimizer(const std::string& filename);

	std::vector<std::string> optimize(int maxTime);

private:
    const Config& config_;

    std::vector<std::string> optimizeTime();
    std::vector<std::string> optimizeStockTargets();

    bool canExecuteProcess(const Process& process, std::vector<Stock>& currentStocks);
    void executeProcess(const Process& process, std::vector<Stock>& currentStocks);
};

 