#include "Parser.hpp"
#include "Stock.hpp"
#include "Process.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <regex>
#include <algorithm>

/**
 * @brief Parse a solution file
 * @param filename The file to parse
 * @return A vector of process executions (cycle, process name)
 */
std::vector<std::pair<int, std::string>> parseSolution(const std::string& filename) {
    std::vector<std::pair<int, std::string>> result;
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return result;
    }
    
    std::string line;
    std::regex lineRegex(R"((\d+):([^:]+))");
    
    while (std::getline(file, line)) {
        std::smatch matches;
        if (std::regex_match(line, matches, lineRegex)) {
            int cycle = std::stoi(matches[1].str());
            std::string processName = matches[2].str();
            
            // Trim whitespace from process name
            processName.erase(0, processName.find_first_not_of(" \t"));
            processName.erase(processName.find_last_not_of(" \t") + 1);
            
            result.emplace_back(cycle, processName);
        }
    }
    
    return result;
}

int main(int argc, char* argv[]) {
    // Check arguments
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <file> <result_to_test>" << std::endl;
        return 1;
    }
    
    std::string configFile = argv[1];
    std::string solutionFile = argv[2];
    
    // Parse the configuration file
    auto parseResult = Parser::parseFile(configFile);
    if (!parseResult) {
        std::cerr << "Error: Failed to parse file " << configFile << std::endl;
        return 1;
    }
    
    const auto& initialStocks = parseResult->stocks;
    const auto& processes = parseResult->processes;
    
    // Create a map of process name to process
    std::unordered_map<std::string, std::shared_ptr<Process>> processMap;
    for (const auto& process : processes) {
        processMap[process->getName()] = process;
    }
    
    // Parse the solution file
    auto solution = parseSolution(solutionFile);
    if (solution.empty()) {
        std::cerr << "Error: Failed to parse solution file " << solutionFile << std::endl;
        return 1;
    }
    
    // Sort solution by cycle
    std::sort(solution.begin(), solution.end(),
              [](const std::pair<int, std::string>& a, const std::pair<int, std::string>& b) {
                  return a.first < b.first;
              });
    
    // Verify the solution
    std::unordered_map<std::string, int> stocks;
    for (const auto& [name, stock] : initialStocks) {
        stocks[name] = stock->getQuantity();
    }
    
    // Map of process name -> end cycle
    std::unordered_map<std::string, int> processEndCycles;
    
    bool valid = true;
    int lastCycle = 0;
    
    for (const auto& [cycle, processName] : solution) {
        lastCycle = std::max(lastCycle, cycle);
        
        // Process any completions
        auto it = processEndCycles.begin();
        while (it != processEndCycles.end()) {
            if (it->second <= cycle) {
                auto processIt = processMap.find(it->first);
                if (processIt != processMap.end()) {
                    // Add results
                    for (const auto& [stockName, quantity] : processIt->second->getResults()) {
                        stocks[stockName] += quantity;
                    }
                }
                it = processEndCycles.erase(it);
            } else {
                ++it;
            }
        }
        
        // Check if the process exists
        auto processIt = processMap.find(processName);
        if (processIt == processMap.end()) {
            std::cerr << "Error at cycle " << cycle << ": Unknown process " << processName << std::endl;
            valid = false;
            continue;
        }
        
        const auto& process = processIt->second;
        
        // Check if we have enough resources
        for (const auto& [stockName, requiredQuantity] : process->getRequirements()) {
            auto stockIt = stocks.find(stockName);
            if (stockIt == stocks.end() || stockIt->second < requiredQuantity) {
                std::cerr << "Error at cycle " << cycle << ": Not enough " << stockName 
                          << " for process " << processName 
                          << " (have " << (stockIt == stocks.end() ? 0 : stockIt->second)
                          << ", need " << requiredQuantity << ")" << std::endl;
                valid = false;
            }
        }
        
        // Consume resources
        for (const auto& [stockName, requiredQuantity] : process->getRequirements()) {
            stocks[stockName] -= requiredQuantity;
        }
        
        // Schedule completion
        processEndCycles[processName] = cycle + process->getCycleDuration();
    }
    
    // Process any remaining completions
    for (const auto& [processName, endCycle] : processEndCycles) {
        lastCycle = std::max(lastCycle, endCycle);
        
        auto processIt = processMap.find(processName);
        if (processIt != processMap.end()) {
            // Add results
            for (const auto& [stockName, quantity] : processIt->second->getResults()) {
                stocks[stockName] += quantity;
            }
        }
    }
    
    // Print result
    if (valid) {
        std::cout << "The solution is valid!" << std::endl;
    } else {
        std::cout << "The solution is INVALID!" << std::endl;
    }
    
    std::cout << "Last cycle: " << lastCycle << std::endl;
    std::cout << "Final stocks:" << std::endl;
    for (const auto& [name, quantity] : stocks) {
        std::cout << "  " << name << " => " << quantity << std::endl;
    }
    
    return valid ? 0 : 1;
}