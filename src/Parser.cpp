#include "Parser.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <regex>
#include <algorithm>

std::optional<Parser::ParseResult> Parser::parseFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return std::nullopt;
    }
    
    ParseResult result;
    std::string line;
    
    while (std::getline(file, line)) {
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') {
            continue;
        }
        
        // Trim leading and trailing whitespace
        line.erase(0, line.find_first_not_of(" \t"));
        line.erase(line.find_last_not_of(" \t") + 1);
        
        if (line.empty()) {
            continue;
        }
        
        // Try to parse as a stock definition
        if (auto stockDef = parseStock(line)) {
            const auto& [name, quantity] = *stockDef;
            result.stocks[name] = std::make_shared<Stock>(name, quantity);
            continue;
        }
        
        // Try to parse as a process definition
        if (auto process = parseProcess(line)) {
            result.processes.push_back(*process);
            continue;
        }
        
        // Try to parse as an optimization target
        if (auto optimization = parseOptimization(line)) {
            result.optimizationTarget = *optimization;
            continue;
        }
        
        std::cerr << "Warning: Ignoring unrecognized line: " << line << std::endl;
    }
    
    return result;
}

std::optional<std::pair<std::string, int>> Parser::parseStock(const std::string& line) {
    std::regex stockRegex(R"(([^:]+):(\d+))");
    std::smatch matches;
    
    if (std::regex_match(line, matches, stockRegex)) {
        std::string name = matches[1].str();
        int quantity = std::stoi(matches[2].str());
        
        // Trim whitespace from name
        name.erase(0, name.find_first_not_of(" \t"));
        name.erase(name.find_last_not_of(" \t") + 1);
        
        return std::make_pair(name, quantity);
    }
    
    return std::nullopt;
}

std::optional<std::shared_ptr<Process>> Parser::parseProcess(const std::string& line) {
    // Format: name:(need1:qty1;need2:qty2;...):(result1:qty1;result2:qty2;...):delay
    std::regex processRegex(R"(([^:]+):\(([^)]*)\):\(([^)]*)\):(\d+))");
    std::smatch matches;
    
    if (std::regex_match(line, matches, processRegex)) {
        std::string name = matches[1].str();
        std::string requirements = matches[2].str();
        std::string results = matches[3].str();
        int delay = std::stoi(matches[4].str());
        
        // Trim whitespace from name
        name.erase(0, name.find_first_not_of(" \t"));
        name.erase(name.find_last_not_of(" \t") + 1);
        
        auto process = std::make_shared<Process>(name, delay);
        
        // Parse requirements
        auto requirementsMap = parseResourceList(requirements);
        for (const auto& [stockName, quantity] : requirementsMap) {
            process->addRequirement(stockName, quantity);
        }
        
        // Parse results
        auto resultsMap = parseResourceList(results);
        for (const auto& [stockName, quantity] : resultsMap) {
            process->addResult(stockName, quantity);
        }
        
        return process;
    }
    
    return std::nullopt;
}

std::optional<Parser::OptimizationTarget> Parser::parseOptimization(const std::string& line) {
    // Format: optimize:(time|stock1;time|stock2;...)
    std::regex optimizeRegex(R"(optimize:\(([^)]*)\))");
    std::smatch matches;
    
    if (std::regex_match(line, matches, optimizeRegex)) {
        std::string targetsList = matches[1].str();
        
        OptimizationTarget target;
        target.optimizeTime = false;
        
        std::istringstream iss(targetsList);
        std::string item;
        
        while (std::getline(iss, item, ';')) {
            // Trim whitespace
            item.erase(0, item.find_first_not_of(" \t"));
            item.erase(item.find_last_not_of(" \t") + 1);
            
            if (item == "time") {
                target.optimizeTime = true;
            } else {
                target.optimizeStocks.push_back(item);
            }
        }
        
        return target;
    }
    
    return std::nullopt;
}

std::unordered_map<std::string, int> Parser::parseResourceList(const std::string& resourceList) {
    std::unordered_map<std::string, int> result;
    
    std::istringstream iss(resourceList);
    std::string item;
    
    while (std::getline(iss, item, ';')) {
        // Skip empty items
        if (item.empty()) {
            continue;
        }
        
        // Parse item in format "name:quantity"
        size_t colonPos = item.find(':');
        if (colonPos != std::string::npos) {
            std::string name = item.substr(0, colonPos);
            int quantity = std::stoi(item.substr(colonPos + 1));
            
            // Trim whitespace from name
            name.erase(0, name.find_first_not_of(" \t"));
            name.erase(name.find_last_not_of(" \t") + 1);
            
            result[name] = quantity;
        }
    }
    
    return result;
}