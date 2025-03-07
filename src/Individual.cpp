#include "Individual.hpp"
#include <limits>

Individual::Individual(const std::vector<int> &sequence)
	: processSequence(sequence), fitness(std::numeric_limits<double>::lowest()) {}