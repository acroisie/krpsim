#pragma once

#include <unordered_map>
#include <vector>
#include <memory>
#include <functional>
#include <string>
#include <chrono>
#include "Solution.hpp"
#include "Parser.hpp"

class Stock;
class Process;
class Simulator;

/**
 * @class Optimizer
 * @brief Finds optimal solutions for process execution
 */
class Optimizer {
public:
    /**
     * @brief Constructor
     * @param stocks The initial stocks
     * @param processes The available processes
     * @param optimizationTarget The optimization targets
     */
    Optimizer(
        const std::unordered_map<std::string, std::shared_ptr<Stock>>& stocks,
        const std::vector<std::shared_ptr<Process>>& processes,
        const Parser::OptimizationTarget& optimizationTarget);
    
    /**
     * @brief Find an optimal solution
     * @param timeLimit The maximum time to spend optimizing (in seconds)
     * @return The best solution found
     */
    Solution findOptimalSolution(int timeLimit);
    
private:
    std::unordered_map<std::string, std::shared_ptr<Stock>> stocks_;
    std::vector<std::shared_ptr<Process>> processes_;
    Parser::OptimizationTarget optimizationTarget_;
    std::unique_ptr<Simulator> simulator_;
    
    /**
     * @brief Evaluate a solution
     * @param solution The solution to evaluate
     * @return A score (higher is better)
     */
    double evaluateSolution(const Solution& solution) const;
    
    /**
     * @brief Generate an initial greedy solution
     * @return A greedy solution
     */
    Solution generateGreedySolution() const;
    
    /**
     * @brief Try to improve a solution using local search
     * @param solution The solution to improve
     * @param timeLimit The time limit for improvement
     * @return The improved solution
     */
    Solution improveWithLocalSearch(Solution solution, const std::chrono::steady_clock::time_point& endTime) const;
};