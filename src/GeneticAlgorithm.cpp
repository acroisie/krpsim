#include "GeneticAlgorithm.hpp"
#include <algorithm>
#include <iostream>
#include <numeric>
#include <unordered_map>

using namespace std;

GeneticAlgorithm::GeneticAlgorithm(const Config &config, Simulator &simulator,
                                   int populationSize, double mutationRate,
                                   double crossoverRate, int eliteCount,
                                   int minSequenceLength, int maxSequenceLength)
    : config(config), simulator(simulator), populationSize(populationSize),
      mutationRate(mutationRate), crossoverRate(crossoverRate),
      eliteCount(min(max(0, eliteCount), populationSize)),
      minSequenceLength(minSequenceLength),
      maxSequenceLength(maxSequenceLength), randomGenerator(random_device{}()) {
    for (const auto &proc : config.getProcesses()) {
        processNames.push_back(proc.name);
    }
}

bool GeneticAlgorithm::canExecuteProcess(const Process* process, const std::map<std::string, int>& stocks) const {
    if (!process) {
        return false;
    }

    for (const auto& [resource, quantity] : process->inputs) {
        auto stockIt = stocks.find(resource);
        if (stockIt == stocks.end() || stockIt->second < quantity) {
            return false;
        }
    }
    return true;
}

void GeneticAlgorithm::updateStocksAfterProcess(const Process* process, std::map<std::string, int>& stocks) const {
    if (!process) {
        return;
    }
    
    // Consume inputs
    for (const auto& [resource, quantity] : process->inputs) {
        stocks[resource] -= quantity;
    }
    
    // Add outputs
    for (const auto& [resource, quantity] : process->outputs) {
        stocks[resource] += quantity;
    }
}

Individual GeneticAlgorithm::createSmartIndividual() {
    // Initialize stocks from config
    map<string, int> currentStocks;
    for (const auto& stock : config.getStocks()) {
        currentStocks[stock.name] = stock.quantity;
    }
    
    // Create a mapping of process names to process pointers for quick lookup
    unordered_map<string, const Process*> processMap;
    for (const auto& process : config.getProcesses()) {
        processMap[process.name] = &process;
    }
    
    // Create sequence
    vector<string> sequence;
    int sequenceLength = 0;
    int maxAttempts = maxSequenceLength * 2; // To avoid potential infinite loops
    int attempts = 0;
    
    while (sequenceLength < maxSequenceLength && attempts < maxAttempts) {
        attempts++;
        
        // Find executable processes
        vector<const Process*> executableProcesses;
        for (const auto& [name, proc] : processMap) {
            if (canExecuteProcess(proc, currentStocks)) {
                executableProcesses.push_back(proc);
            }
        }
        
        if (executableProcesses.empty()) {
            // No executable processes found, break the loop
            break;
        }
        
        // Select a random executable process with preference for processes leading to optimization goals
        uniform_int_distribution<size_t> dist(0, executableProcesses.size() - 1);
        const Process* selectedProcess = executableProcesses[dist(randomGenerator)];
        sequence.push_back(selectedProcess->name);
        sequenceLength++;
        
        // Update stocks
        updateStocksAfterProcess(selectedProcess, currentStocks);
    }
    
    // If sequence is too short, pad it with random processes
    // This helps maintain genetic diversity
    if (sequence.size() < minSequenceLength) {
        uniform_int_distribution<size_t> procDist(0, processNames.size() - 1);
        while (sequence.size() < minSequenceLength) {
            sequence.push_back(processNames[procDist(randomGenerator)]);
        }
    }
    
    return Individual(sequence);
}

Individual GeneticAlgorithm::createRandomIndividual() {
    uniform_int_distribution<> lengthDist(minSequenceLength, maxSequenceLength);
    uniform_int_distribution<> processDist(
        0, static_cast<int>(processNames.size()) - 1);

    int sequenceLength = lengthDist(randomGenerator);
    vector<string> sequence;
    sequence.reserve(sequenceLength);

    for (int i = 0; i < sequenceLength; ++i) {
        sequence.push_back(processNames[processDist(randomGenerator)]);
    }

    return Individual(sequence);
}

void GeneticAlgorithm::initializePopulation() {
    population.clear();
    population.reserve(populationSize);

    // Create a mix of smart and random individuals
    // Smart individuals will guide the search, random ones maintain diversity
    int smartCount = populationSize * 0.8; // 80% smart, 20% random
    
    for (int i = 0; i < smartCount; ++i) {
        population.push_back(createSmartIndividual());
    }
    
    for (int i = smartCount; i < populationSize; ++i) {
        population.push_back(createRandomIndividual());
    }

    cout << "Population initialized with " << populationSize << " individuals: " 
         << smartCount << " smart, " << (populationSize - smartCount) << " random." << endl;
}

void GeneticAlgorithm::evaluatePopulation() {
    for (auto &individual : population) {
        Simulator::SimulationResult result =
            simulator.runSimulation(individual.processSequence);
        individual.fitnessScore = result.fitness;
    }
}

vector<size_t> GeneticAlgorithm::selectParents() {
    vector<size_t> selectedIndices;
    if (population.empty()) {
        return selectedIndices;
    }

    // Find valid min and max fitness values
    double minFitness = numeric_limits<double>::max();
    double maxFitness = numeric_limits<double>::lowest();

    for (const auto &individual : population) {
        if (isfinite(individual.fitnessScore)) {
            minFitness = min(minFitness, individual.fitnessScore);
            maxFitness = max(maxFitness, individual.fitnessScore);
        }
    }

    // Handle edge cases
    if (minFitness == numeric_limits<double>::max() ||
        maxFitness == numeric_limits<double>::lowest()) {
        // No valid fitness scores, select randomly
        selectedIndices.reserve(populationSize);
        uniform_int_distribution<size_t> randDist(0, population.size() - 1);
        for (int i = 0; i < populationSize; ++i) {
            selectedIndices.push_back(randDist(randomGenerator));
        }
        return selectedIndices;
    }

    // If all fitness values are identical
    if (minFitness == maxFitness) {
        // Select parents randomly
        selectedIndices.reserve(populationSize);
        uniform_int_distribution<size_t> randDist(0, population.size() - 1);
        for (int i = 0; i < populationSize; ++i) {
            selectedIndices.push_back(randDist(randomGenerator));
        }
        return selectedIndices;
    }

    // Normalize fitness scores to positive values
    double range = maxFitness - minFitness;
    vector<double> normalizedFitness;
    normalizedFitness.reserve(population.size());
    double totalNormalizedFitness = 0.0;

    for (const auto &individual : population) {
        double normalized = 0.001; // Minimum score
        if (isfinite(individual.fitnessScore)) {
            normalized = (individual.fitnessScore - minFitness) / range + 0.001;
        }
        normalizedFitness.push_back(normalized);
        totalNormalizedFitness += normalized;
    }

    // Handle invalid total fitness
    if (totalNormalizedFitness <= 0.0 || !isfinite(totalNormalizedFitness)) {
        // Select parents randomly
        selectedIndices.reserve(populationSize);
        uniform_int_distribution<size_t> randDist(0, population.size() - 1);
        for (int i = 0; i < populationSize; ++i) {
            selectedIndices.push_back(randDist(randomGenerator));
        }
        return selectedIndices;
    }

    // Roulette wheel selection
    uniform_real_distribution<> dist(0.0, totalNormalizedFitness);
    selectedIndices.reserve(populationSize);

    for (int i = 0; i < populationSize; ++i) {
        double randomPoint = dist(randomGenerator);
        double currentSum = 0.0;
        bool selected = false;

        for (size_t j = 0; j < population.size(); ++j) {
            currentSum += normalizedFitness[j];
            if (randomPoint <= currentSum) {
                selectedIndices.push_back(j);
                selected = true;
                break;
            }
        }

        // Fallback if no selection made
        if (!selected && !population.empty()) {
            selectedIndices.push_back(population.size() - 1);
        }
    }

    return selectedIndices;
}

pair<Individual, Individual>
GeneticAlgorithm::crossover(const Individual &parent1,
                            const Individual &parent2) {
    uniform_real_distribution<> crossDist(0.0, 1.0);
    if (crossDist(randomGenerator) > crossoverRate) {
        return {parent1, parent2};
    }

    const auto &seq1 = parent1.processSequence;
    const auto &seq2 = parent2.processSequence;
    size_t len1 = seq1.size();
    size_t len2 = seq2.size();

    if (len1 < 2 || len2 < 2) {
        return {parent1, parent2};
    }

    // Two-point crossover
    size_t shortLength = min(len1, len2);
    uniform_int_distribution<size_t> pointDist(0, shortLength - 1);

    size_t point1 = pointDist(randomGenerator);
    size_t point2 = pointDist(randomGenerator);

    if (point1 == point2) {
        point2 = (point1 + 1) % shortLength;
    }

    if (point1 > point2) {
        swap(point1, point2);
    }

    vector<string> childSeq1, childSeq2;
    childSeq1.reserve(len1);
    childSeq2.reserve(len2);

    // Create first child
    childSeq1.insert(childSeq1.end(), seq1.begin(), seq1.begin() + point1);
    if (point2 < len2) {
        childSeq1.insert(childSeq1.end(), seq2.begin() + point1,
                         seq2.begin() + point2);
    } else if (point1 < len2) {
        childSeq1.insert(childSeq1.end(), seq2.begin() + point1, seq2.end());
    }
    if (point2 < len1) {
        childSeq1.insert(childSeq1.end(), seq1.begin() + point2, seq1.end());
    }

    // Create second child
    childSeq2.insert(childSeq2.end(), seq2.begin(), seq2.begin() + point1);
    if (point2 < len1) {
        childSeq2.insert(childSeq2.end(), seq1.begin() + point1,
                         seq1.begin() + point2);
    } else if (point1 < len1) {
        childSeq2.insert(childSeq2.end(), seq1.begin() + point1, seq1.end());
    }
    if (point2 < len2) {
        childSeq2.insert(childSeq2.end(), seq2.begin() + point2, seq2.end());
    }

    return {Individual(childSeq1), Individual(childSeq2)};
}

Individual GeneticAlgorithm::mutate(const Individual &individual) {
    Individual mutated = individual;
    uniform_real_distribution<> probDist(0.0, 1.0);
    uniform_int_distribution<> processDist(
        0, static_cast<int>(processNames.size()) - 1);

    // Point mutations
    for (size_t i = 0; i < mutated.processSequence.size(); ++i) {
        if (probDist(randomGenerator) < mutationRate) {
            mutated.processSequence[i] =
                processNames[processDist(randomGenerator)];
        }
    }

    // Structural mutations (add/remove elements)
    double structuralMutationRate = 0.02;
    if (probDist(randomGenerator) < structuralMutationRate &&
        mutated.processSequence.size() > 1) {
        uniform_int_distribution<size_t> posDist(
            0, mutated.processSequence.size());
        size_t pos = posDist(randomGenerator);

        // 50% chance to remove or add an element
        if (probDist(randomGenerator) < 0.5 &&
            mutated.processSequence.size() >
                static_cast<size_t>(minSequenceLength)) {
            // Remove element if we have more than the minimum
            if (pos < mutated.processSequence.size()) {
                mutated.processSequence.erase(mutated.processSequence.begin() +
                                              pos);
            }
        } else {
            // Add element
            string randomProcess = processNames[processDist(randomGenerator)];
            if (pos > mutated.processSequence.size()) {
                pos = mutated.processSequence.size();
            }
            mutated.processSequence.insert(
                mutated.processSequence.begin() + pos, randomProcess);
        }
    }

    return mutated;
}

void GeneticAlgorithm::selectNextGeneration() {
    vector<Individual> newPopulation;
    newPopulation.reserve(populationSize);

    // Sort population by fitness (descending)
    sort(population.begin(), population.end());

    // Keep elite individuals
    for (int i = 0; i < eliteCount && i < static_cast<int>(population.size());
         ++i) {
        newPopulation.push_back(population[i]);
    }

    // Select parents for breeding
    vector<size_t> parentIndices = selectParents();
    if (parentIndices.empty() && populationSize > eliteCount) {
        // If no parents selected, fill with elites or random individuals
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

    // Create next generation
    uniform_int_distribution<size_t> parentDist(0, parentIndices.size() - 1);
    while (newPopulation.size() < static_cast<size_t>(populationSize)) {
        size_t idx1 = parentDist(randomGenerator);
        size_t idx2 = parentDist(randomGenerator);

        if (idx1 >= parentIndices.size() || idx2 >= parentIndices.size()) {
            continue;
        }

        size_t parentIndex1 = parentIndices[idx1];
        size_t parentIndex2 = parentIndices[idx2];

        if (parentIndex1 >= population.size() ||
            parentIndex2 >= population.size()) {
            continue;
        }

        // Try to avoid self-crossover
        int tries = 0;
        while (parentIndex1 == parentIndex2 && parentIndices.size() > 1 &&
               tries < 10) {
            idx2 = parentDist(randomGenerator);
            if (idx2 >= parentIndices.size()) continue;
            parentIndex2 = parentIndices[idx2];
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

    for (int gen = 0; gen < generations; ++gen) {
        evaluatePopulation();
        Individual currentBest = getBestIndividual();

        // Update overall best
        if (gen == 0 || currentBest.fitnessScore > bestOverall.fitnessScore) {
            bestOverall = currentBest;
            cout << "Generation " << (gen + 1) << "/" << generations
                 << " - New Best Fitness: " << bestOverall.fitnessScore << endl;
        }

        selectNextGeneration();
    }

    // Final evaluation
    evaluatePopulation();
    Individual finalBest = getBestIndividual();
    if (finalBest.fitnessScore > bestOverall.fitnessScore) {
        bestOverall = finalBest;
    }

    cout << "Evolution finished. Final best fitness: "
         << bestOverall.fitnessScore << endl;
    return bestOverall;
}