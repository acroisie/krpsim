#pragma once
#include <string>
#include <vector>
#include <map>
using namespace std;

struct Process
{
	string name;
	map<string, string> inputs;
	map<string, string> outputs;
	unsigned int time;
};

class Parser
{
public:
	Parser(string &filename);
	~Parser();

	const map<string, int> &getStocks() const;
	const vector<Process> &getProcesses() const;
	const vector<string> &getOptimizeGoals() const;

private:
	void parseLine(string &line);

	map<string, int> stocks;
	vector<Process> processes;
	vector<string> optimizeGoals;
	enum class Section { STOCKS, PROCESSES, OPTIMIZE_GOALS };
	Section currentSection;
};