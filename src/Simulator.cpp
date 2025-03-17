#include "Simulator.hpp"
#include "Stock.hpp"
#include "Process.hpp"
#include <algorithm>
#include <iostream>
#include <set>

Simulator::Simulator(
    const std::unordered_map<std::string, std::shared_ptr<Stock>>& stocks,
    const std::vector<std::shared_ptr<Process>>& processes)
    : initialStocks_(stocks), processes_(processes) {
    
    // Create process map for quick lookup
    for (const auto& process : processes_) {
        processMap_[process->getName()] = process;
    }
}

std::unordered_map<int, Simulator::StockState> Simulator::simulate(const Solution& solution) const {
    std::unordered_map<int, StockState> result;
    
    // Initialize with initial stocks
    StockState initialState;
    for (const auto& [name, stock] : initialStocks_) {
        initialState.quantities[name] = stock->getQuantity();
    }
    result[0] = initialState;
    
    // Get all execution cycles (start and end)
    std::set<int> cycles;
    cycles.insert(0);
    for (const auto& execution : solution.getExecutions()) {
        cycles.insert(execution.startCycle);
        cycles.insert(execution.endCycle);
    }
    
    // Simulate cycle by cycle
    int prevCycle = 0;
    for (int cycle : cycles) {
        if (cycle == 0) continue;
        
        // Copy previous state
        result[cycle] = result[prevCycle];
        
        // Process completions
        for (const auto& execution : solution.getExecutions()) {
            if (execution.endCycle == cycle) {
                auto processIt = processMap_.find(execution.processName);
                if (processIt != processMap_.end()) {
                    // Add results
                    for (const auto& [stockName, quantity] : processIt->second->getResults()) {
                        result[cycle].quantities[stockName] += quantity;
                    }
                }
            }
        }
        
        // Process starts
        for (const auto& execution : solution.getExecutions()) {
            if (execution.startCycle == cycle) {
                auto processIt = processMap_.find(execution.processName);
                if (processIt != processMap_.end()) {
                    // Consume requirements
                    for (const auto& [stockName, quantity] : processIt->second->getRequirements()) {
                        result[cycle].quantities[stockName] -= quantity;
                    }
                }
            }
        }
        
        prevCycle = cycle;
    }
    
    return result;
}

bool Simulator::validate(const Solution& solution) const {
    auto cycleStates = simulate(solution);
    
    // Check that no stock goes negative
    for (const auto& [cycle, state] : cycleStates) {
        for (const auto& [stockName, quantity] : state.quantities) {
            if (quantity < 0) {
                std::cerr << "Error: Negative stock " << stockName 
                          << " with quantity " << quantity
                          << " at cycle " << cycle << std::endl;
                return false;
            }
        }
    }
    
    return true;
}

std::vector<std::string> Simulator::findExecutableProcesses(const StockState& stockState) const {
    std::vector<std::string> result;
    
    // Create a map of stocks
    std::unordered_map<std::string, std::shared_ptr<Stock>> tempStocks;
    for (const auto& [name, quantity] : stockState.quantities) {
        tempStocks[name] = std::make_shared<Stock>(name, quantity);
    }
    
    // Check each process
    for (const auto& process : processes_) {
        if (process->canExecute(tempStocks)) {
            result.push_back(process->getName());
        }
    }
    
    return result;
}