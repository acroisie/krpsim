#include "Individual.hpp"
#include <limits>

Individual::Individual() {}
Individual::Individual(const std::vector<std::string> &sequence)
	: processSequence(sequence), fitness(std::numeric_limits<double>::lowest()) {}