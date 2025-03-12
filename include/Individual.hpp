#pragma once
#include <limits>
#include <string>
#include <vector>

class Individual {
  public:
    Individual();
    explicit Individual(const std::vector<std::string> &processSequence);

    std::vector<std::string> processSequence;
    double fitness;
};
