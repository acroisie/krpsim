#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include "Stock.hpp"
#include "Process.hpp"

class Parser
{
public:
	Parser(const std::string &filename);
	bool parse();

	const std::vector<Stock> &getStocks() const;
	const std::vector<Process> &getProcesses() const;

private:
	std::string _filename;
	std::vector<Stock> _stocks;
	std::vector<Process> _processes;

	void parseStock(const std::string &line);
	void parseProcess(const std::string &line);
	void parseResources(const std::string &token, std::unordered_map<std::string, int> &resources);
	void parseOptimization(const std::string &line);
};