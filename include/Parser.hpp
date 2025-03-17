#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <memory>
#include <optional>
#include "Process.hpp"
#include "Stock.hpp"

/**
 * @class Parser
 * @brief Parses configuration files for the krpsim program
 */
class Parser {
public:
    /**
     * @struct OptimizationTarget
     * @brief Represents what should be optimized
     */
    struct OptimizationTarget {
        bool optimizeTime;
        std::vector<std::string> optimizeStocks;
    };
    
    /**
     * @struct ParseResult
     * @brief Contains the result of parsing a configuration file
     */
    struct ParseResult {
        std::unordered_map<std::string, std::shared_ptr<Stock>> stocks;
        std::vector<std::shared_ptr<Process>> processes;
        OptimizationTarget optimizationTarget;
    };
    
    /**
     * @brief Parse a configuration file
     * @param filename The file to parse
     * @return ParseResult containing stocks, processes, and optimization targets, or std::nullopt on error
     */
    static std::optional<ParseResult> parseFile(const std::string& filename);
    
private:
    static std::optional<std::pair<std::string, int>> parseStock(const std::string& line);
    static std::optional<std::shared_ptr<Process>> parseProcess(const std::string& line);
    static std::optional<OptimizationTarget> parseOptimization(const std::string& line);
    static std::unordered_map<std::string, int> parseResourceList(const std::string& resourceList);
};