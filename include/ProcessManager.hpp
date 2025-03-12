#pragma once
#include "Config.hpp"
#include "Individual.hpp"
#include <map>
#include <string>
#include <vector>

class ProcessManager {
  public:
    ProcessManager(const Config &config, int delayLimit);
    void runGeneticAlgorithm();

  private:
    const Config &config_;
    int delayLimit_;
    int currentCycle_;

    std::map<std::string, int> currentStocks_;
    std::vector<std::pair<int, std::string>> executionLogs_;

    std::vector<const Process *> getRunnableProcesses();
    bool executeProcess(const Process *process);

    void initializePopulation();
    void evaluateFitness();

    double calculateFitness(Individual &individual);
    std::vector<Individual> selectParents(const std::vector<Individual> &population);
    std::pair<Individual, Individual> crossover(const Individual &parent1, const Individual &parent2);
    Individual mutate(const Individual &individual, double mutationRate);

    Individual findBestIndividual(const std::vector<Individual> &population);

    void generateOutput();

    static const int POPULATION_SIZE;
    std::vector<Individual> population_;
    Individual bestSolution_;
};
