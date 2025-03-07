#pragma once
#include <vector>
#include <string>

class Individual {
public:
	Individual(const std::vector<int> &processSequence);

	std::vector<int> processSequence;
	double fitness;
};
