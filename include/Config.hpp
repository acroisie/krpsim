#pragma once
#include <string>
#include <vector>
#include "Stock.hpp"
#include "Process.hpp"

class Config {
public:
	explicit Config(const std::string &filename);
	bool parse();

	const std::vector<Stock> &getStocks() const;
	const std::vector<Process> &getProcesses() const;
	const std::vector<std::string> &getOptimizationTargets() const;

	bool isValid() const;

private:
	std::vector<Stock> stocks_;
	std::vector<Process> processes_;
	std::vector<std::string> optimizationTargets_;

	bool isValid;

	void parseConfigFile(const std::string &filename);
};