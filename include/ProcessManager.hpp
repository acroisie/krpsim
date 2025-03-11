#pragma once
#include <vector>
#include <string>
#include <map>
#include "Config.hpp"
#include "Individual.hpp"

class ProcessManager {
public:
    ProcessManager(const Config &config, int delayLimit);
    bool runGeneticAlgorithm();

private:
    const Config &config_;
    int delayLimit_;
    int currentCycle_;

    std::map<std::string, int> currentStocks_;
    std::vector<std::pair<int, std::string>> executionLogs_;

    std::vector<const Process*> getRunnableProcesses();
    bool executeProcess(const Process* process);
    void generateOutput();
    Individual mutate(const Individual& individual, double mutationRate);
    std::vector<Individual> mutationPopulation(std::vector<Individual>& population, double mutationRate);
    void initializePopulation();
    void evaluateFitness();
    std::pair<Individual, Individual> crossover(const Individual& parent1, const Individual& parent2);
    std::vector<Individual> crossoverPopulation(const std::vector<Individual>& parents, double crossoverRate);
    std::vector<Individual> selectParents(const std::vector<Individual>& population);


    static const int POPULATION_SIZE;
    std::vector<Individual> population_;
    std::vector<Individual> selection();

    double calculateFitness(Individual &individual);
};
