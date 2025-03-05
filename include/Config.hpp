#pragma once
#include <string>
#include <vector>
#include "Stock.hpp"
#include "Process.hpp"

class Config{
public:
	void addStock(const Stock &stock);
	void addProcess(const Process &process);
	void setOptimizeGoal(const std::vector<std::string> &goal);

	const std::vector<Stock> &getStocks() const;
	const std::vector<Process> &getProcesses() const;
	const std::vector<std::string> &getOptimizeGoal() const;

private:
	std::vector<Stock> stocks_;
	std::vector<Process> processes_;
	std::vector<std::string> optimizeGoal_;
};