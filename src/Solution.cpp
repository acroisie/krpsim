#include "Solution.hpp"
#include "Stock.hpp"
#include "Process.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>

Solution::Solution() {}

void Solution::addExecution(const std::string& processName, int startCycle, int duration) {
    executions_.emplace_back(processName, startCycle, startCycle + duration);
}

const std::vector<Solution::ProcessExecution>& Solution::getExecutions() const {
    return executions_;
}

std::vector<std::string> Solution::getProcessesStartingAt(int cycle) const {
    std::vector<std::string> result;
    
    for (const auto& execution : executions_) {
        if (execution.startCycle == cycle) {
            result.push_back(execution.processName);
        }
    }
    
    return result;
}

int Solution::getTotalDuration() const {
    if (executions_.empty()) {
        return 0;
    }
    
    return std::max_element(
        executions_.begin(),
        executions_.end(),
        [](const ProcessExecution& a, const ProcessExecution& b) {
            return a.endCycle < b.endCycle;
        }
    )->endCycle;
}

std::unordered_map<std::string, int> Solution::getFinalStocks(
    const std::unordered_map<std::string, std::shared_ptr<Stock>>& initialStocks,
    const std::vector<std::shared_ptr<Process>>& processes) const {
    
    // Create a map of process name to process
    std::unordered_map<std::string, std::shared_ptr<Process>> processMap;
    for (const auto& process : processes) {
        processMap[process->getName()] = process;
    }
    
    // Create a copy of the initial stocks
    std::unordered_map<std::string, int> finalStocks;
    for (const auto& [name, stock] : initialStocks) {
        finalStocks[name] = stock->getQuantity();
    }
    
    // Sort executions by start cycle
    std::vector<ProcessExecution> sortedExecutions = executions_;
    std::sort(sortedExecutions.begin(), sortedExecutions.end(),
              [](const ProcessExecution& a, const ProcessExecution& b) {
                  return a.startCycle < b.startCycle;
              });
    
    // Simulate the solution
    for (const auto& execution : sortedExecutions) {
        auto processIt = processMap.find(execution.processName);
        if (processIt == processMap.end()) {
            std::cerr << "Error: Unknown process " << execution.processName << std::endl;
            continue;
        }
        
        const auto& process = processIt->second;
        
        // Consume requirements
        for (const auto& [stockName, quantity] : process->getRequirements()) {
            finalStocks[stockName] -= quantity;
            if (finalStocks[stockName] < 0) {
                std::cerr << "Error: Negative stock " << stockName << " at cycle " 
                          << execution.startCycle << " for process " << execution.processName << std::endl;
            }
        }
        
        // Produce results
        for (const auto& [stockName, quantity] : process->getResults()) {
            finalStocks[stockName] += quantity;
        }
    }
    
    return finalStocks;
}

bool Solution::saveToFile(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing" << std::endl;
        return false;
    }
    
    // Sort executions by start cycle
    std::vector<ProcessExecution> sortedExecutions = executions_;
    std::sort(sortedExecutions.begin(), sortedExecutions.end(),
              [](const ProcessExecution& a, const ProcessExecution& b) {
                  return a.startCycle < b.startCycle;
              });
    
    for (const auto& execution : sortedExecutions) {
        file << execution.startCycle << ":" << execution.processName << "\n";
    }
    
    return true;
}