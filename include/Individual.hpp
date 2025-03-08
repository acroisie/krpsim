#pragma once
#include <vector>
#include <string>

class Individual {
public:
	Individual();
	Individual(const std::vector<std::string> &processSequence);

	std::vector<std::string> processSequence;
	double fitness;
};
