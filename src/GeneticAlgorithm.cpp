#include "GeneticAlgorithm.hpp"
#include <algorithm>
#include <iostream>
#include <numeric>
#include <unordered_map>

using namespace std;

GeneticAlgorithm::GeneticAlgorithm(const Config &configRef,
                                   Simulator &simulatorRef, int populationSize,
                                   double mutationRate, double crossoverRate,
                                   int eliteCount, int minSequenceLength,
                                   int maxSequenceLength)
    : config(configRef), simulator(simulatorRef),
      populationSize(populationSize), mutationRate(mutationRate),
      crossoverRate(crossoverRate),
      eliteCount(min(max(0, eliteCount), populationSize)),
      minSequenceLength(minSequenceLength),
      maxSequenceLength(maxSequenceLength), randomGenerator(random_device{}()) {
    for (const auto &process : config.getProcesses())
        processNames.push_back(process.name);
}

bool GeneticAlgorithm::canExecuteProcess(
    const Process *process, const std::map<std::string, int> &stocks) const {
    if (!process) {
        return false;
    }

    for (const auto &[resource, quantity] : process->inputs) {
        auto stockIt = stocks.find(resource);
        if (stockIt == stocks.end() || stockIt->second < quantity) {
            return false;
        }
    }
    return true;
}

void GeneticAlgorithm::updateStocksAfterProcess(
    const Process *process, std::map<std::string, int> &stocks) const {
    if (!process) {
        return;
    }

    for (const auto &[resource, quantity] : process->inputs) {
        stocks[resource] -= quantity;
    }

    for (const auto &[resource, quantity] : process->outputs) {
        stocks[resource] += quantity;
    }
}

Individual GeneticAlgorithm::createSmartIndividual() {
    map<string, int> availableStocks;
    for (const auto &stock : config.getStocks())
        availableStocks[stock.name] = stock.quantity;

    unordered_map<string, const Process *> processMap;
    for (const auto &process : config.getProcesses())
        processMap[process.name] = &process;

    vector<string> processSequence;
    int attempts = 0;
    while (processSequence.size() < (static_cast<size_t>(maxSequenceLength)) &&
           attempts < maxSequenceLength * 2) {
        ++attempts;
        vector<const Process *> executableProcesses;
        for (const auto &[name, processPtr] : processMap)
            if (canExecuteProcess(processPtr, availableStocks))
                executableProcesses.push_back(processPtr);

        if (executableProcesses.empty()) break;

        auto bestProcess = *min_element(
            executableProcesses.begin(), executableProcesses.end(),
            [&](const Process *processA, const Process *processB) {
                const auto &priorityMap = simulator.getProcessPriority();
                int priorityA = priorityMap.count(processA->name)
                                    ? priorityMap.at(processA->name)
                                    : 3;
                int priorityB = priorityMap.count(processB->name)
                                    ? priorityMap.at(processB->name)
                                    : 3;
                if (priorityA != priorityB) return priorityA < priorityB;
                return processA->nbCycle < processB->nbCycle;
            });
        processSequence.push_back(bestProcess->name);
        updateStocksAfterProcess(bestProcess, availableStocks);
    }

    if (processSequence.size() < static_cast<size_t>(minSequenceLength) &&
        !processSequence.empty()) {
        std::mt19937 &rng = randomGenerator;
        while (processSequence.size() < static_cast<size_t>(minSequenceLength))
            processSequence.push_back(
                processSequence[rng() % processSequence.size()]);
    }

    return Individual(processSequence);
}

Individual GeneticAlgorithm::createRandomIndividual() {
    uniform_int_distribution<> lengthDistribution(minSequenceLength,
                                                  maxSequenceLength);
    uniform_int_distribution<> processIndexDistribution(
        0, static_cast<int>(processNames.size()) - 1);

    int processSequenceLength = lengthDistribution(randomGenerator);
    vector<string> processSequence;
    processSequence.reserve(processSequenceLength);

    for (int processIndex = 0; processIndex < processSequenceLength;
         ++processIndex) {
        processSequence.push_back(
            processNames[processIndexDistribution(randomGenerator)]);
    }

    return Individual(processSequence);
}

void GeneticAlgorithm::initializePopulation() {
    population.clear();
    population.reserve(populationSize);

    int smartIndividualCount = populationSize * 0.8;

    for (int i = 0; i < smartIndividualCount; ++i) {
        population.push_back(createSmartIndividual());
    }

    for (int i = smartIndividualCount; i < populationSize; ++i) {
        population.push_back(createRandomIndividual());
    }

    cout << "Population initialized with " << populationSize
         << " individuals: " << smartIndividualCount << " smart, "
         << (populationSize - smartIndividualCount) << " random." << endl;
}

void GeneticAlgorithm::evaluatePopulation() {
    for (auto &individual : population) {
        Simulator::SimulationResult simulationResult =
            simulator.runSimulation(individual.processSequence);
        individual.fitnessScore = simulationResult.fitness;
    }
}

vector<size_t> GeneticAlgorithm::selectParents() {
    vector<size_t> selectedParentIndices;
    if (population.empty()) {
        return selectedParentIndices;
    }

    double minFitness = numeric_limits<double>::max();
    double maxFitness = numeric_limits<double>::lowest();

    for (const auto &individual : population) {
        if (isfinite(individual.fitnessScore)) {
            minFitness = min(minFitness, individual.fitnessScore);
            maxFitness = max(maxFitness, individual.fitnessScore);
        }
    }

    if (minFitness == numeric_limits<double>::max() ||
        maxFitness == numeric_limits<double>::lowest()) {
        selectedParentIndices.reserve(populationSize);
        uniform_int_distribution<size_t> randomIndexDistribution(
            0, population.size() - 1);
        for (int i = 0; i < populationSize; ++i) {
            selectedParentIndices.push_back(
                randomIndexDistribution(randomGenerator));
        }
        return selectedParentIndices;
    }

    if (minFitness == maxFitness) {
        selectedParentIndices.reserve(populationSize);
        uniform_int_distribution<size_t> randomIndexDistribution(
            0, population.size() - 1);
        for (int i = 0; i < populationSize; ++i) {
            selectedParentIndices.push_back(
                randomIndexDistribution(randomGenerator));
        }
        return selectedParentIndices;
    }

    double fitnessRange = maxFitness - minFitness;
    vector<double> normalizedFitness;
    normalizedFitness.reserve(population.size());
    double totalNormalizedFitness = 0.0;

    for (const auto &individual : population) {
        double normalized = 0.001;
        if (isfinite(individual.fitnessScore)) {
            normalized =
                (individual.fitnessScore - minFitness) / fitnessRange + 0.001;
        }
        normalizedFitness.push_back(normalized);
        totalNormalizedFitness += normalized;
    }

    if (totalNormalizedFitness <= 0.0 || !isfinite(totalNormalizedFitness)) {
        selectedParentIndices.reserve(populationSize);
        uniform_int_distribution<size_t> randomIndexDistribution(
            0, population.size() - 1);
        for (int i = 0; i < populationSize; ++i) {
            selectedParentIndices.push_back(
                randomIndexDistribution(randomGenerator));
        }
        return selectedParentIndices;
    }

    uniform_real_distribution<> randomFitnessDistribution(
        0.0, totalNormalizedFitness);
    selectedParentIndices.reserve(populationSize);

    for (int i = 0; i < populationSize; ++i) {
        double randomPoint = randomFitnessDistribution(randomGenerator);
        double currentSum = 0.0;
        bool parentSelected = false;

        for (size_t index = 0; index < population.size(); ++index) {
            currentSum += normalizedFitness[index];
            if (randomPoint <= currentSum) {
                selectedParentIndices.push_back(index);
                parentSelected = true;
                break;
            }
        }

        if (!parentSelected && !population.empty()) {
            selectedParentIndices.push_back(population.size() - 1);
        }
    }

    return selectedParentIndices;
}

pair<Individual, Individual>
GeneticAlgorithm::crossover(const Individual &parent1,
                            const Individual &parent2) {
    uniform_real_distribution<> crossoverProbabilityDistribution(0.0, 1.0);
    if (crossoverProbabilityDistribution(randomGenerator) > crossoverRate) {
        return {parent1, parent2};
    }

    const auto &sequence1 = parent1.processSequence;
    const auto &sequence2 = parent2.processSequence;
    size_t length1 = sequence1.size();
    size_t length2 = sequence2.size();

    if (length1 < 2 || length2 < 2) {
        return {parent1, parent2};
    }

    size_t shortestLength = min(length1, length2);
    uniform_int_distribution<size_t> crossoverPointDistribution(
        0, shortestLength - 1);

    size_t crossoverPoint1 = crossoverPointDistribution(randomGenerator);
    size_t crossoverPoint2 = crossoverPointDistribution(randomGenerator);

    if (crossoverPoint1 == crossoverPoint2) {
        crossoverPoint2 = (crossoverPoint1 + 1) % shortestLength;
    }

    if (crossoverPoint1 > crossoverPoint2) {
        swap(crossoverPoint1, crossoverPoint2);
    }

    vector<string> childSequence1, childSequence2;
    childSequence1.reserve(length1);
    childSequence2.reserve(length2);

    childSequence1.insert(childSequence1.end(), sequence1.begin(),
                          sequence1.begin() + crossoverPoint1);
    if (crossoverPoint2 < length2) {
        childSequence1.insert(childSequence1.end(),
                              sequence2.begin() + crossoverPoint1,
                              sequence2.begin() + crossoverPoint2);
    } else if (crossoverPoint1 < length2) {
        childSequence1.insert(childSequence1.end(),
                              sequence2.begin() + crossoverPoint1,
                              sequence2.end());
    }
    if (crossoverPoint2 < length1) {
        childSequence1.insert(childSequence1.end(),
                              sequence1.begin() + crossoverPoint2,
                              sequence1.end());
    }

    childSequence2.insert(childSequence2.end(), sequence2.begin(),
                          sequence2.begin() + crossoverPoint1);
    if (crossoverPoint2 < length1) {
        childSequence2.insert(childSequence2.end(),
                              sequence1.begin() + crossoverPoint1,
                              sequence1.begin() + crossoverPoint2);
    } else if (crossoverPoint1 < length1) {
        childSequence2.insert(childSequence2.end(),
                              sequence1.begin() + crossoverPoint1,
                              sequence1.end());
    }
    if (crossoverPoint2 < length2) {
        childSequence2.insert(childSequence2.end(),
                              sequence2.begin() + crossoverPoint2,
                              sequence2.end());
    }

    return {Individual(childSequence1), Individual(childSequence2)};
}

Individual GeneticAlgorithm::mutate(const Individual &individual) {
    Individual mutatedIndividual = individual;
    uniform_real_distribution<> mutationProbabilityDistribution(0.0, 1.0);
    uniform_int_distribution<> processIndexDistribution(
        0, static_cast<int>(processNames.size()) - 1);

    for (size_t processIndex = 0;
         processIndex < mutatedIndividual.processSequence.size();
         ++processIndex) {
        if (mutationProbabilityDistribution(randomGenerator) < mutationRate) {
            mutatedIndividual.processSequence[processIndex] =
                processNames[processIndexDistribution(randomGenerator)];
        }
    }

    double structuralMutationRate = 0.02;
    if (mutationProbabilityDistribution(randomGenerator) <
            structuralMutationRate &&
        mutatedIndividual.processSequence.size() > 1) {
        uniform_int_distribution<size_t> positionDistribution(
            0, mutatedIndividual.processSequence.size());
        size_t position = positionDistribution(randomGenerator);

        if (mutationProbabilityDistribution(randomGenerator) < 0.5 &&
            mutatedIndividual.processSequence.size() >
                static_cast<size_t>(minSequenceLength)) {
            if (position < mutatedIndividual.processSequence.size()) {
                mutatedIndividual.processSequence.erase(
                    mutatedIndividual.processSequence.begin() + position);
            }
        } else {
            string randomProcess =
                processNames[processIndexDistribution(randomGenerator)];
            if (position > mutatedIndividual.processSequence.size()) {
                position = mutatedIndividual.processSequence.size();
            }
            mutatedIndividual.processSequence.insert(
                mutatedIndividual.processSequence.begin() + position,
                randomProcess);
        }
    }

    return mutatedIndividual;
}

void GeneticAlgorithm::selectNextGeneration() {
    vector<Individual> newPopulation;
    newPopulation.reserve(populationSize);

    sort(population.begin(), population.end());

    for (int i = 0; i < eliteCount && i < static_cast<int>(population.size());
         ++i) {
        newPopulation.push_back(population[i]);
    }

    vector<size_t> parentIndices = selectParents();
    if (parentIndices.empty() && populationSize > eliteCount) {
        while (newPopulation.size() < static_cast<size_t>(populationSize)) {
            if (!population.empty()) {
                newPopulation.push_back(population[0]);
            } else {
                newPopulation.push_back(createSmartIndividual());
            }
        }
        population = newPopulation;
        return;
    }

    if (parentIndices.empty()) {
        population = newPopulation;
        return;
    }

    uniform_int_distribution<size_t> parentIndexDistribution(
        0, parentIndices.size() - 1);
    while (newPopulation.size() < static_cast<size_t>(populationSize)) {
        size_t index1 = parentIndexDistribution(randomGenerator);
        size_t index2 = parentIndexDistribution(randomGenerator);

        if (index1 >= parentIndices.size() || index2 >= parentIndices.size()) {
            continue;
        }

        size_t parentIndex1 = parentIndices[index1];
        size_t parentIndex2 = parentIndices[index2];

        if (parentIndex1 >= population.size() ||
            parentIndex2 >= population.size()) {
            continue;
        }

        int tries = 0;
        while (parentIndex1 == parentIndex2 && parentIndices.size() > 1 &&
               tries < 10) {
            index2 = parentIndexDistribution(randomGenerator);
            if (index2 >= parentIndices.size()) continue;
            parentIndex2 = parentIndices[index2];
            if (parentIndex2 >= population.size()) continue;
            tries++;
        }

        const Individual &parent1 = population[parentIndex1];
        const Individual &parent2 = population[parentIndex2];

        auto [child1, child2] = crossover(parent1, parent2);
        child1 = mutate(child1);
        child2 = mutate(child2);

        if (newPopulation.size() < static_cast<size_t>(populationSize)) {
            newPopulation.push_back(child1);
        }
        if (newPopulation.size() < static_cast<size_t>(populationSize)) {
            newPopulation.push_back(child2);
        }
    }

    population = newPopulation;
}

Individual GeneticAlgorithm::getBestIndividual() const {
    if (population.empty()) {
        return Individual();
    }

    auto bestIt = max_element(population.begin(), population.end(),
                              [](const Individual &a, const Individual &b) {
                                  return a.fitnessScore < b.fitnessScore;
                              });
    return *bestIt;
}

Individual GeneticAlgorithm::runEvolution(int generations) {
    cout << "Starting genetic algorithm evolution for " << generations
         << " generations..." << endl;

    initializePopulation();
    Individual bestOverall;
    bestOverall.fitnessScore = numeric_limits<double>::lowest();

    for (int generation = 0; generation < generations; ++generation) {
        evaluatePopulation();
        Individual currentBest = getBestIndividual();

        if (generation == 0 ||
            currentBest.fitnessScore > bestOverall.fitnessScore) {
            bestOverall = currentBest;
            cout << "Generation " << (generation + 1) << "/" << generations
                 << " - New Best Fitness: " << bestOverall.fitnessScore << endl;
        }

        selectNextGeneration();
    }

    evaluatePopulation();
    Individual finalBest = getBestIndividual();
    if (finalBest.fitnessScore > bestOverall.fitnessScore) {
        bestOverall = finalBest;
    }

    cout << "Evolution finished. Final best fitness: "
         << bestOverall.fitnessScore << endl;
    return bestOverall;
}