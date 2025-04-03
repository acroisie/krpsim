#pragma once

#include <limits>
#include <string>
#include <vector>

class Individual {
  public:
    std::vector<std::string> processSequence;
    double fitnessScore;

    Individual() : fitnessScore(std::numeric_limits<double>::lowest()) {}

    bool operator<(const Individual &other) const {
        return fitnessScore < other.fitnessScore;
    }

    bool operator>(const Individual &other) const {
        return fitnessScore > other.fitnessScore;
    }
};