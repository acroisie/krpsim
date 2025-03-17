#pragma once

#include <string>
#include <unordered_map>
#include <memory>
#include <vector>

class Stock;

/**
 * @class Process
 * @brief Represents a process that consumes and produces stocks
 */
class Process {
public:
    /**
     * @brief Constructor
     * @param name The name of the process
     * @param cycleDuration The number of cycles this process takes
     */
    Process(const std::string& name, int cycleDuration);
    
    /**
     * @brief Add an input requirement for this process
     * @param stockName The name of the required stock
     * @param quantity The quantity needed
     */
    void addRequirement(const std::string& stockName, int quantity);
    
    /**
     * @brief Add an output result for this process
     * @param stockName The name of the produced stock
     * @param quantity The quantity produced
     */
    void addResult(const std::string& stockName, int quantity);
    
    /**
     * @brief Get the name of the process
     * @return The process name
     */
    const std::string& getName() const;
    
    /**
     * @brief Get the cycle duration of the process
     * @return The cycle duration
     */
    int getCycleDuration() const;
    
    /**
     * @brief Check if the process can be executed with the given stocks
     * @param stocks Map of stock name to stock pointer
     * @return True if the process can be executed
     */
    bool canExecute(const std::unordered_map<std::string, std::shared_ptr<Stock>>& stocks) const;
    
    /**
     * @brief Execute the process on the given stocks
     * @param stocks Map of stock name to stock pointer
     * @return True if successful
     */
    bool execute(std::unordered_map<std::string, std::shared_ptr<Stock>>& stocks);
    
    /**
     * @brief Get the requirements of this process
     * @return Map of stock name to required quantity
     */
    const std::unordered_map<std::string, int>& getRequirements() const;
    
    /**
     * @brief Get the results of this process
     * @return Map of stock name to produced quantity
     */
    const std::unordered_map<std::string, int>& getResults() const;
    
private:
    std::string name_;
    int cycleDuration_;
    std::unordered_map<std::string, int> requirements_;
    std::unordered_map<std::string, int> results_;
};