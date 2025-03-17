#pragma once

#include <unordered_map>
#include <vector>
#include <memory>
#include <string>
#include "Solution.hpp"

class Stock;
class Process;

/**
 * @class Simulator
 * @brief Simulates the execution of processes
 */
class Simulator {
public:
    /**
     * @struct StockState
     * @brief Represents the state of stocks at a specific cycle
     */
    struct StockState {
        std::unordered_map<std::string, int> quantities;
    };
    
    /**
     * @brief Constructor
     * @param stocks The initial stocks
     * @param processes The available processes
     */
    Simulator(
        const std::unordered_map<std::string, std::shared_ptr<Stock>>& stocks,
        const std::vector<std::shared_ptr<Process>>& processes);
    
    /**
     * @brief Simulate a solution
     * @param solution The solution to simulate
     * @return Map of cycle to stock state
     */
    std::unordered_map<int, StockState> simulate(const Solution& solution) const;
    
    /**
     * @brief Simulate and validate a solution
     * @param solution The solution to validate
     * @return True if the solution is valid
     */
    bool validate(const Solution& solution) const;
    
    /**
     * @brief Find processes that can be executed at a given cycle
     * @param stockState The current stock state
     * @return Vector of executable process names
     */
    std::vector<std::string> findExecutableProcesses(const StockState& stockState) const;
    
private:
    std::unordered_map<std::string, std::shared_ptr<Stock>> initialStocks_;
    std::vector<std::shared_ptr<Process>> processes_;
    std::unordered_map<std::string, std::shared_ptr<Process>> processMap_;
};