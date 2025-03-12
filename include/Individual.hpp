#pragma once
#include <vector>
#include <string>
#include <limits>

class Individual {
public:
    Individual();
    explicit Individual(const std::vector<std::string> &processSequence);

    std::vector<std::string> processSequence;
    double fitness;
};
