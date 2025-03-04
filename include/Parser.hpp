#pragma once
#include <string>
#include <vector>
#include "Stock.hpp"
#include "Process.hpp"

class Parser {
public:
    explicit Parser(const std::string& filename);
    bool parse();

    const std::vector<Stock>& getStocks() const;
    const std::vector<Process>& getProcesses() const;

private:
    std::string filename_;
    std::vector<Stock> stocks_;
    std::vector<Process> processes_;
    bool optimizeFound_;

    void parseLine(const std::string& line, int lineNumber);
};
