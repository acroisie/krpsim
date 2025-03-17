#include "Optimizer.hpp"
#include "Stock.hpp"
#include "Process.hpp"
#include "Simulator.hpp"
#include <algorithm>
#include <iostream>
#include <random>
#include <set>

Optimizer::Optimizer(
    const std::unordered_map<std::string, std::shared_ptr<Stock>>& stocks,
    const std::vector<std::shared_ptr<Process>>& processes,
    const Parser::OptimizationTarget& optimizationTarget)
    : stocks_(stocks), processes_(processes), optimizationTarget_(optimizationTarget) {
    
    simulator_ = std::make_unique<Simulator>(stocks, processes);
}

Solution Optimizer::findOptimalSolution(int timeLimit) {
    auto startTime = std::chrono::steady_clock::now();
    auto endTime = startTime + std::chrono::seconds(timeLimit);
    
    std::cout << "Finding optimal solution (time limit: " << timeLimit << " seconds)..." << std::endl;
    
    // Generate an initial greedy solution
    Solution bestSolution = generateGreedySolution();
    double bestScore = evaluateSolution(bestSolution);
    
    std::cout << "Initial solution score: " << bestScore << std::endl;
    
    // Improve the solution using local search
    Solution improvedSolution = improveWithLocalSearch(bestSolution, endTime);
    double improvedScore = evaluateSolution(improvedSolution);
    
    if (improvedScore > bestScore) {
        bestSolution = improvedSolution;
        bestScore = improvedScore;
        std::cout << "Improved solution score: " << bestScore << std::endl;
    }
    
    auto elapsed = std::chrono::steady_clock::now() - startTime;
    std::cout << "Optimization completed in " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count() / 1000.0
              << " seconds" << std::endl;
    
    return bestSolution;
}

double Optimizer::evaluateSolution(const Solution& solution) const {
    // Validate the solution
    if (!simulator_->validate(solution)) {
        return -1000000.0;  // Invalid solution gets a very bad score
    }
    
    double score = 0.0;
    
    // Get final stocks
    auto finalStocks = solution.getFinalStocks(stocks_, processes_);
    
    if (optimizationTarget_.optimizeTime) {
        // Lower duration is better
        score -= solution.getTotalDuration();
    }
    
    // Maximize target stocks
    for (const auto& stockName : optimizationTarget_.optimizeStocks) {
        auto it = finalStocks.find(stockName);
        if (it != finalStocks.end()) {
            score += it->second;
        }
    }
    
    return score;
}

Solution Optimizer::generateGreedySolution() const {
    Solution solution;
    
    // Copy initial stocks
    std::unordered_map<std::string, int> currentStocks;
    for (const auto& [name, stock] : stocks_) {
        currentStocks[name] = stock->getQuantity();
    }
    
    // Keep track of the current cycle
    int currentCycle = 0;
    
    // Map process name to its end cycle
    std::unordered_map<std::string, int> processEndCycles;
    
    // Set of processes that can be executed at the current cycle
    std::set<std::string> executableProcesses;
    
    // Create temporary stocks for checking
    std::unordered_map<std::string, std::shared_ptr<Stock>> tempStocks;
    for (const auto& [name, quantity] : currentStocks) {
        tempStocks[name] = std::make_shared<Stock>(name, quantity);
    }
    
    // Create a process map for quick lookup
    std::unordered_map<std::string, std::shared_ptr<Process>> processMap;
    for (const auto& process : processes_) {
        processMap[process->getName()] = process;
    }
    
    // Find initially executable processes
    for (const auto& process : processes_) {
        if (process->canExecute(tempStocks)) {
            executableProcesses.insert(process->getName());
        }
    }
    
    // Maximum number of iterations to prevent infinite loops
    const int MAX_ITERATIONS = 10000;
    int iterations = 0;
    
    while (!executableProcesses.empty() && iterations < MAX_ITERATIONS) {
        iterations++;
        
        // Choose a process to execute
        std::string processToExecute = *executableProcesses.begin();
        
        // For optimization, prioritize processes that produce target stocks
        for (const auto& processName : executableProcesses) {
            const auto& process = processMap[processName];
            
            // Check if this process produces any target stocks
            bool producesTargetStock = false;
            for (const auto& targetStock : optimizationTarget_.optimizeStocks) {
                if (process->getResults().find(targetStock) != process->getResults().end()) {
                    producesTargetStock = true;
                    break;
                }
            }
            
            if (producesTargetStock) {
                processToExecute = processName;
                break;
            }
        }
        
        // Execute the process
        const auto& process = processMap[processToExecute];
        
        // Consume requirements
        for (const auto& [stockName, quantity] : process->getRequirements()) {
            currentStocks[stockName] -= quantity;
        }
        
        // Add to solution
        int duration = process->getCycleDuration();
        solution.addExecution(processToExecute, currentCycle, duration);
        
        // Update end cycle
        processEndCycles[processToExecute] = currentCycle + duration;
        
        // Update stocks when the process completes
        for (const auto& [stockName, quantity] : process->getResults()) {
            // The stock will be available after the process completes
            int endCycle = processEndCycles[processToExecute];
            
            // We'll add the result when we reach that cycle
            if (currentCycle == endCycle) {
                currentStocks[stockName] += quantity;
            }
        }
        
        // Remove executed process from the set
        executableProcesses.erase(processToExecute);
        
        // Update executable processes
        tempStocks.clear();
        for (const auto& [name, quantity] : currentStocks) {
            tempStocks[name] = std::make_shared<Stock>(name, quantity);
        }
        
        for (const auto& process : processes_) {
            if (process->canExecute(tempStocks) && 
                executableProcesses.find(process->getName()) == executableProcesses.end()) {
                executableProcesses.insert(process->getName());
            }
        }
        
        // If no executable processes, advance to the next process completion
        if (executableProcesses.empty() && !processEndCycles.empty()) {
            int nextCycle = std::numeric_limits<int>::max();
            for (const auto& [process, endCycle] : processEndCycles) {
                if (endCycle > currentCycle && endCycle < nextCycle) {
                    nextCycle = endCycle;
                }
            }
            
            if (nextCycle != std::numeric_limits<int>::max()) {
                currentCycle = nextCycle;
                
                // Update stocks for completed processes
                for (auto it = processEndCycles.begin(); it != processEndCycles.end(); ) {
                    if (it->second == currentCycle) {
                        // Process completed, add results
                        const auto& process = processMap[it->first];
                        for (const auto& [stockName, quantity] : process->getResults()) {
                            currentStocks[stockName] += quantity;
                        }
                        
                        // Remove from end cycles
                        it = processEndCycles.erase(it);
                    } else {
                        ++it;
                    }
                }
                
                // Update executable processes
                tempStocks.clear();
                for (const auto& [name, quantity] : currentStocks) {
                    tempStocks[name] = std::make_shared<Stock>(name, quantity);
                }
                
                for (const auto& process : processes_) {
                    if (process->canExecute(tempStocks)) {
                        executableProcesses.insert(process->getName());
                    }
                }
            }
        }
    }
    
    return solution;
}

Solution Optimizer::improveWithLocalSearch(Solution solution, const std::chrono::steady_clock::time_point& endTime) const {
    // This is a simple local search implementation that could be expanded
    
    Solution bestSolution = solution;
    double bestScore = evaluateSolution(bestSolution);
    
    // Random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    
    while (std::chrono::steady_clock::now() < endTime) {
        // Make a small modification to the current solution
        Solution newSolution = bestSolution;
        
        // Get the executions
        auto executions = newSolution.getExecutions();
        
        if (executions.empty()) {
            break;
        }
        
        // Choose a random execution to modify
        std::uniform_int_distribution<> execDist(0, executions.size() - 1);
        // int execIndex = execDist(gen);  // Uncomment when implementing modification logic
        
        // Currently we just move the execution earlier or later in time
        std::uniform_int_distribution<> timeDist(-10, 10);
        // int timeShift = timeDist(gen);  // Uncomment when implementing modification logic
        
        // TODO: Implement logic to modify solution based on execIndex and timeShift
        
        // Modify the solution (this is a simplified approach)
        // For a real implementation, we would need to rebuild the solution properly
        
        // Evaluate the new solution
        double newScore = evaluateSolution(newSolution);
        
        // If the new solution is better, keep it
        if (newScore > bestScore) {
            bestSolution = newSolution;
            bestScore = newScore;
        }
    }
    
    return bestSolution;
}