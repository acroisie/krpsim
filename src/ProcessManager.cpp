#include "ProcessManager.hpp"
#include <iostream>
#include <limits>
#include <vector>
#include <string>
#include <algorithm>
#include <random>

using namespace std;

const int ProcessManager::POPULATION_SIZE = 50;

ProcessManager::ProcessManager(const Config &config, int delayLimit)
    : config_(config), delayLimit_(delayLimit), currentCycle_(0) {
    for (const auto &stock : config_.getStocks()) {
        currentStocks_[stock.name] = stock.quantity;
    }
}

bool ProcessManager::runSimulation() {
    while (currentCycle_ < delayLimit_) {
        vector<const Process*> runnableProcesses = getRunnableProcesses();
        if (runnableProcesses.empty()) {
            break;
        }

        currentCycle_++;
    }
    generateOutput();
    return true;
}

void ProcessManager::initializePopulation() {
    population_.resize(POPULATION_SIZE);
    vector<string> availableProcessNames;
    for (const auto &process : config_.getProcesses()) {
        availableProcessNames.push_back(process.name);
    }

    random_device rd;
    mt19937 generator(rd());
    uniform_int_distribution<> lengthDistribution(1, delayLimit_);
    uniform_int_distribution<> processDistribution(0, availableProcessNames.size() - 1);

    for (int i = 0; i < POPULATION_SIZE; ++i) {
        vector<string> randomSequence;
        int sequenceLength = lengthDistribution(generator);
        for (int j = 0; j < sequenceLength; ++j) {
            if (!availableProcessNames.empty()) {
                int processIndex = processDistribution(generator);
                randomSequence.push_back(availableProcessNames[processIndex]);
            }
        }
        population_[i] = Individual(randomSequence);
    }
}

void ProcessManager::evaluateFitness() {
    for (int i = 0; i < POPULATION_SIZE; ++i) {
        population_[i].fitness = calculateFitness(population_[i]);
    }
}

vector<Individual> ProcessManager::selection() {
    vector<Individual> parents = selectParents(population_);
    cout << "Parents selected." << endl;
    return parents;
}

bool ProcessManager::runGeneticAlgorithm() {
    double mutationRate = 0.1;
    double crossoverRate = 0.7;
    int generationCount = 50;

    cout << "Running genetic algorithm... (runGeneticAlgorithm() START)" << endl;

    initializePopulation();
    cout << "Initial population generated." << endl;

    cout << "Starting evolution over " << generationCount << " generations..." << endl;
    for (int generation = 0; generation < generationCount; ++generation) {
        cout << "--- Generation " << generation << " ---" << endl;

        evaluateFitness();
        vector<Individual> parents = selection();
        vector<Individual> nextGeneration = crossoverPopulation(parents, crossoverRate);
        population_ = mutationPopulation(nextGeneration, mutationRate);

        double bestFitnessThisGen = numeric_limits<double>::lowest();
        for (const auto& indiv : population_) {
            if (indiv.fitness > bestFitnessThisGen) {
                bestFitnessThisGen = indiv.fitness;
            }
        }
        cout << "  Best fitness in generation " << generation << ": " << bestFitnessThisGen << endl;
    }
    cout << "Evolution finished over " << generationCount << " generations." << endl;

    evaluateFitness();
    generateOutput();

    cout << "Running genetic algorithm... (runGeneticAlgorithm() END)" << endl;
    return true;
}

double ProcessManager::calculateFitness(Individual &individual) {
    currentStocks_.clear();
    for (const auto &stock : config_.getStocks()) {
        currentStocks_[stock.name] = stock.quantity;
    }
    executionLogs_.clear();
    currentCycle_ = 0;

    for (const string& processName : individual.processSequence) {
        if (currentCycle_ >= delayLimit_) {
            break;
        }
        bool processExecutedThisCycle = false;

        do {
            processExecutedThisCycle = false;
            vector<const Process*> runnableProcesses = getRunnableProcesses();

            if (!runnableProcesses.empty()) {
                const Process* processToExecute = nullptr;
                for (const auto &proc : runnableProcesses) {
                    if (proc->name == processName) {
                        processToExecute = proc;
                        break;
                    }
                }
                if (!processToExecute && !runnableProcesses.empty()) {
                    processToExecute = runnableProcesses[0];
                }

                if (processToExecute) {
                    executeProcess(processToExecute);
                    executionLogs_.push_back({currentCycle_, processToExecute->name});
                    processExecutedThisCycle = true;
                } else {
                    break;
                }
            } else {
                break;
            }

        } while (processExecutedThisCycle);
        currentCycle_++;
    }

    double fitnessScore = 0.0;
    const vector<string>& optimizeGoals = config_.getOptimizeGoal();
    if (!optimizeGoals.empty()) {
        string goal = optimizeGoals[0];
        if (currentStocks_.count(goal) > 0) {
            fitnessScore = currentStocks_[goal];
        } else if (goal == "time") {
            fitnessScore = -currentCycle_;
        } else {
            fitnessScore = numeric_limits<double>::lowest();
        }
    } else {
        fitnessScore = 0.0;
    }
    return fitnessScore;
}

vector<Individual> ProcessManager::selectParents(const vector<Individual>& population) {
    vector<Individual> parents;
    vector<double> fitnessValues;
    for (const auto& individual : population) {
        fitnessValues.push_back(individual.fitness);
    }

    double totalFitness = 0.0;
    for (double fitness : fitnessValues) {
        totalFitness += fitness;
    }

    vector<double> selectionProbabilities;
    for (double fitness : fitnessValues) {
        selectionProbabilities.push_back(fitness / totalFitness);
    }

    random_device rd;
    mt19937 generator(rd());
    uniform_real_distribution<> distribution(0.0, 1.0);

    for (int i = 0; i < POPULATION_SIZE; ++i) { // Sélectionne POPULATION_SIZE parents pour la prochaine génération
        double randomValue = distribution(generator);
        double cumulativeProbability = 0.0;
        for (size_t j = 0; j < selectionProbabilities.size(); ++j) {
            cumulativeProbability += selectionProbabilities[j];
            if (randomValue <= cumulativeProbability) {
                parents.push_back(population[j]);
                // break;
            }
        }
    }
    return parents;
}

pair<Individual, Individual> ProcessManager::crossover(const Individual& parent1, const Individual& parent2) {
    random_device rd;
    mt19937 generator(rd());
    uniform_int_distribution<> distribution(0, min(parent1.processSequence.size(), parent2.processSequence.size()) > 0 ? min(parent1.processSequence.size(), parent2.processSequence.size()) - 1 : 0); // Prevent distribution from being invalid if sequences are empty

    size_t crossoverPoint = 0;
    if (parent1.processSequence.size() > 0 && parent2.processSequence.size() > 0) { // Only crossover if sequences are not empty
        crossoverPoint = distribution(generator);
    }

    vector<string> child1Sequence;
    vector<string> child2Sequence;

    for (size_t i = 0; i < parent1.processSequence.size(); ++i) {
        if (i <= crossoverPoint) {
            child1Sequence.push_back(parent1.processSequence[i]);
        } else if (i < parent2.processSequence.size()){
            child1Sequence.push_back(parent2.processSequence[i]);
        }
    }

    for (size_t i = 0; i < parent2.processSequence.size(); ++i) {
         if (i <= crossoverPoint && i < parent1.processSequence.size()) {
            child2Sequence.push_back(parent2.processSequence[i]);
        } else if (i < parent1.processSequence.size()){
            child2Sequence.push_back(parent1.processSequence[i]);
        }
    }

    return make_pair(Individual(child1Sequence), Individual(child2Sequence));
}

vector<Individual> ProcessManager::crossoverPopulation(const vector<Individual>& parents, double crossoverRate) {
    vector<Individual> nextGeneration;
    random_device rd;
    mt19937 generator(rd());
    uniform_real_distribution<> probabilityDistribution(0.0, 1.0);

    for (int i = 0; i < POPULATION_SIZE / 2; ++i) {
        if (probabilityDistribution(generator) < crossoverRate) {
            pair<Individual, Individual> children = crossover(parents[i % parents.size()], parents[(i + 1) % parents.size()]);
            nextGeneration.push_back(children.first);
            nextGeneration.push_back(children.second);
        } else {
            nextGeneration.push_back(parents[i % parents.size()]);
            nextGeneration.push_back(parents[(i + 1) % parents.size()]);
        }
    }
    return nextGeneration;
}

Individual ProcessManager::mutate(const Individual& individual, double mutationRate) {
    random_device rd;
    mt19937 generator(rd());
    uniform_real_distribution<> probabilityDistribution(0.0, 1.0);
    uniform_int_distribution<> processDistribution(0, config_.getProcesses().size() - 1);

    vector<string> mutatedSequence = individual.processSequence;

    if (probabilityDistribution(generator) < mutationRate) {
        if (!mutatedSequence.empty()) {
            uniform_int_distribution<> pointDistribution(0, mutatedSequence.size() - 1);
            int mutationPoint = pointDistribution(generator);
            int processIndex = processDistribution(generator);
            string processNameToMutateTo = config_.getProcesses()[processIndex].name;
            mutatedSequence[mutationPoint] = processNameToMutateTo;
        }
    }

    return Individual(mutatedSequence);
}

vector<Individual> ProcessManager::mutationPopulation(vector<Individual>& population, double mutationRate) {
    vector<Individual> mutatedPopulation;
    for (Individual& individual : population) {
        mutatedPopulation.push_back(mutate(individual, mutationRate));
    }
    return mutatedPopulation;
}

vector<const Process*> ProcessManager::getRunnableProcesses() {
    vector<const Process*> runnable;
    for (const auto &process : config_.getProcesses()) {
        bool canRun = true;
        for (const auto &input : process.inputs) {
            if (currentStocks_.count(input.first) == 0) {
                canRun = false;
                break;
            }
            if (currentStocks_[input.first] < input.second) {
                canRun = false;
                break;
            }
        }
        if (canRun) {
            runnable.push_back(&process);
        }
    }
    return runnable;
}

bool ProcessManager::executeProcess(const Process* process) {
    if (!process) {
        return false;
    }

    for (const auto &input : process->inputs) {
        currentStocks_[input.first] -= input.second;
        if (currentStocks_[input.first] < 0) {
            return false;
        }
    }

    for (const auto &output : process->outputs) {
        currentStocks_[output.first] += output.second;
    }

    return true;
}

void ProcessManager::generateOutput() {
    cout << "Simulation completed in " << currentCycle_ << " cycles." << endl;
    cout << "Final stocks:" << endl;
    for (auto stock : currentStocks_) {
        cout << "Stock " << stock.first << ": " << stock.second << endl;
    }
    cout << "Execution logs:" << endl;
    for (auto log : executionLogs_) {
        cout << "Cycle " << log.first << ": " << log.second << endl;
    }
    // Find and Output Best Individual from the final population
    double bestFitness = numeric_limits<double>::lowest();
    const Individual* bestIndividual = nullptr;
    for (const auto& indiv : population_) {
        if (indiv.fitness > bestFitness) {
            bestFitness = indiv.fitness;
            bestIndividual = &indiv;
        }
    }

    if (bestIndividual) {
        cout << "\nBest Individual in Final Population:" << endl;
        cout << "  Fitness: " << bestIndividual->fitness << endl;
        cout << "  Process Sequence: [";
        for (const auto& processName : bestIndividual->processSequence) {
            cout << processName << ", ";
        }
        cout << "]" << endl;
    } else {
        cout << "\nNo best individual found in final population." << endl;
    }
}