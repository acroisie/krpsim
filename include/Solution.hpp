#pragma once

#include <vector>
#include <string>
#include <unordered_map>
#include <memory>

class Stock;
class Process;

/**
 * @class Solution
 * @brief Represents a solution to the optimization problem
 */
class Solution {
public:
    /**
     * @struct ProcessExecution
     * @brief Represents a single process execution in the schedule
     */
    struct ProcessExecution {
        std::string processName;
        int startCycle;
        int endCycle;
        
        ProcessExecution(const std::string& name, int start, int end)
            : processName(name), startCycle(start), endCycle(end) {}
    };
    
    /**
     * @brief Constructor
     */
    Solution();
    
    /**
     * @brief Add a process execution to the solution
     * @param processName The name of the process
     * @param startCycle The cycle when the process starts
     * @param duration The duration of the process
     */
    void addExecution(const std::string& processName, int startCycle, int duration);
    
    /**
     * @brief Get all process executions in the solution
     * @return Vector of all process executions
     */
    const std::vector<ProcessExecution>& getExecutions() const;
    
    /**
     * @brief Get process executions starting at a specific cycle
     * @param cycle The cycle number
     * @return Vector of process executions starting at the given cycle
     */
    std::vector<std::string> getProcessesStartingAt(int cycle) const;
    
    /**
     * @brief Get the total duration of the solution
     * @return The last cycle in the solution
     */
    int getTotalDuration() const;
    
    /**
     * @brief Get the final stocks after executing the solution
     * @param initialStocks The initial stocks
     * @param processes The available processes
     * @return The final stocks after executing the solution
     */
    std::unordered_map<std::string, int> getFinalStocks(
        const std::unordered_map<std::string, std::shared_ptr<Stock>>& initialStocks,
        const std::vector<std::shared_ptr<Process>>& processes) const;
    
    /**
     * @brief Save the solution to a file
     * @param filename The file to save to
     * @return True if successful
     */
    bool saveToFile(const std::string& filename) const;
    
private:
    std::vector<ProcessExecution> executions_;
};