#include "Parser.hpp"
#include "Optimizer.hpp"
#include "Simulator.hpp"
#include <iostream>
#include <string>
#include <chrono>
#include <iomanip>

int main(int argc, char* argv[]) {
    // Check arguments
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <file> <delay>" << std::endl;
        return 1;
    }
    
    std::string filename = argv[1];
    int timeLimit;
    
    try {
        timeLimit = std::stoi(argv[2]);
    } catch (const std::exception& e) {
        std::cerr << "Error: Invalid delay value" << std::endl;
        return 1;
    }
    
    if (timeLimit <= 0) {
        std::cerr << "Error: Delay must be positive" << std::endl;
        return 1;
    }
    
    // Parse the configuration file
    auto parseResult = Parser::parseFile(filename);
    if (!parseResult) {
        std::cerr << "Error: Failed to parse file " << filename << std::endl;
        return 1;
    }
    
    const auto& stocks = parseResult->stocks;
    const auto& processes = parseResult->processes;
    const auto& optimizationTarget = parseResult->optimizationTarget;
    
    // Print summary
    std::cout << "Nice file! " << processes.size() << " processes, " 
              << stocks.size() << " stocks, " 
              << (optimizationTarget.optimizeTime ? 1 : 0) + optimizationTarget.optimizeStocks.size() 
              << " to optimize" << std::endl;
    
    // Start the optimization
    std::cout << "Evaluating ..................";
    
    Optimizer optimizer(stocks, processes, optimizationTarget);
    Solution solution = optimizer.findOptimalSolution(timeLimit);
    
    std::cout << " done." << std::endl;
    
    // Print the solution
    std::cout << "Main walk" << std::endl;
    
    // Sort executions by start cycle
    std::vector<Solution::ProcessExecution> sortedExecutions = solution.getExecutions();
    std::sort(sortedExecutions.begin(), sortedExecutions.end(),
              [](const Solution::ProcessExecution& a, const Solution::ProcessExecution& b) {
                  return a.startCycle < b.startCycle;
              });
    
    for (const auto& execution : sortedExecutions) {
        std::cout << execution.startCycle << ":" << execution.processName << std::endl;
    }
    
    // Check if we can continue or if we're done
    bool canContinue = false;
    
    auto finalStocks = solution.getFinalStocks(stocks, processes);
    Simulator simulator(stocks, processes);
    Simulator::StockState finalState;
    finalState.quantities = finalStocks;
    
    auto executableProcesses = simulator.findExecutableProcesses(finalState);
    if (executableProcesses.empty()) {
        int finalCycle = solution.getTotalDuration();
        std::cout << "no more process doable at time " << finalCycle << std::endl;
    } else {
        std::cout << "system is self-sustained, stopping at cycle " 
                  << solution.getTotalDuration() << std::endl;
    }
    
    // Print final stocks
    std::cout << "Stock :" << std::endl;
    for (const auto& [name, quantity] : finalStocks) {
        std::cout << "  " << name << " => " << quantity << std::endl;
    }
    
    // Save the solution to a file
    std::string solutionFile = filename + ".solution";
    if (solution.saveToFile(solutionFile)) {
        std::cout << "Solution saved to " << solutionFile << std::endl;
    }
    
    return 0;
}