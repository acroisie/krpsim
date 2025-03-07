#pragma once
#include <vector>
#include <string>

class Individual {
public:
	Individual();
	Individual(const std::vector<int> &processSequence);

	std::vector<int> processSequence;
	double fitness;
};
