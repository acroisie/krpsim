#pragma once

#include <string>
#include <unordered_map>

struct Process
{
	std::string name;
	std::unordered_map<std::string, int> inputs;
	std::unordered_map<std::string, int> outputs;
	int cycleAmount;
};