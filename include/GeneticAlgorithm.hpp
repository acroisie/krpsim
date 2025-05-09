#pragma once
#include "Config.hpp"
#include "Individual.hpp"
#include "Simulator.hpp"
#include <map>
#include <random>
#include <string>
#include <vector>

class GeneticAlgorithm {
  public:
    GeneticAlgorithm(const Config &config, Simulator &simulator,
                     int populationSize = 100, double mutationRate = 0.05,
                     double crossoverRate = 0.7, int eliteCount = 2,
                     int minSequenceLength = 50, int maxSequenceLength = 150);

    Individual runEvolution(int generations);
    Individual getBestIndividual() const;

  private:
    const Config &config;
    Simulator &simulator;
    int populationSize;
    double mutationRate;
    double crossoverRate;
    int eliteCount;
    int minSequenceLength;
    int maxSequenceLength;

    std::vector<Individual> population;
    std::vector<std::string> processNames;
    std::mt19937 randomGenerator;

    void initializePopulation();
    void evaluatePopulation();
    void selectNextGeneration();

    Individual createRandomIndividual();
    Individual createSmartIndividual();

    bool canExecuteProcess(const Process *process,
                           const std::map<std::string, int> &stocks) const;
    void updateStocksAfterProcess(const Process *process,
                                  std::map<std::string, int> &stocks) const;

    std::vector<size_t> selectParents();
    std::pair<Individual, Individual> crossover(const Individual &parent1,
                                                const Individual &parent2);
    Individual mutate(const Individual &individual);
};