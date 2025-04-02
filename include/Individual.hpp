#pragma once
#include <limits>
#include <string>
#include <vector>

class Individual {
  public:
    Individual() : fitnessScore(std::numeric_limits<double>::lowest()) {}

    explicit Individual(const std::vector<std::string> &sequence)
        : processSequence(sequence),
          fitnessScore(std::numeric_limits<double>::lowest()) {}

    std::vector<std::string> processSequence;
    double fitnessScore;

    bool operator<(const Individual &other) const {
        return fitnessScore > other.fitnessScore;
    }

    bool operator>(const Individual &other) const {
        return fitnessScore < other.fitnessScore;
    }
};