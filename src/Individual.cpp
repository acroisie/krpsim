#include "Individual.hpp"

Individual::Individual() : fitness(std::numeric_limits<double>::lowest()) {}

Individual::Individual(const std::vector<std::string> &sequence)
    : processSequence(sequence), fitness(std::numeric_limits<double>::lowest()) {}
