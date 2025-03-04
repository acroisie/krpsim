#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include "Stock.hpp"
#include "Process.hpp"

class Parser
{
public:
	explicit Parser(const std::string &filename);
	bool parse();

	const std::vector<Stock> &getStocks() const;
	const std::vector<Process> &getProcesses() const;

private:
	std::string _filename;
	std::vector<Stock> _stocks;
	std::vector<Process> _processes;

	void parseLine(const std::string &line, int lineNumber);
};