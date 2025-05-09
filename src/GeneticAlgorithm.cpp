#include "GeneticAlgorithm.hpp"
#include <algorithm>
#include <iostream>
#include <numeric>
#include <unordered_map>

using namespace std;

GeneticAlgorithm::GeneticAlgorithm(const Config &cfg, Simulator &sim, int popSize,
                                   double mutRate, double crossRate, int elite,
                                   int minLen, int maxLen)
    : config(cfg), simulator(sim), populationSize(popSize), mutationRate(mutRate),
      crossoverRate(crossRate), eliteCount(min(max(0, elite), popSize)),
      minSequenceLength(minLen), maxSequenceLength(maxLen), randomGenerator(random_device{}()) {
    for (const auto &p : config.getProcesses()) processNames.push_back(p.name);
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
    map<string, int> stocks;
    for (const auto &s : config.getStocks()) stocks[s.name] = s.quantity;

    unordered_map<string, const Process *> pMap;
    for (const auto &p : config.getProcesses()) pMap[p.name] = &p;

    vector<string> seq;
    int attempts = 0;
    while (seq.size() <(static_cast<size_t>(maxSequenceLength)) && attempts < maxSequenceLength * 2) {
        ++attempts;
        vector<const Process *> executable;
        for (const auto &[name, ptr] : pMap)
            if (canExecuteProcess(ptr, stocks)) executable.push_back(ptr);

        if (executable.empty()) break;

        auto best = *min_element(executable.begin(), executable.end(), [&](const Process *a, const Process *b) {
            const auto &prio = simulator.getProcessPriority();
            int pa = prio.count(a->name) ? prio.at(a->name) : 3;
            int pb = prio.count(b->name) ? prio.at(b->name) : 3;
            if (pa != pb) return pa < pb;
            return a->nbCycle < b->nbCycle;
        });
        seq.push_back(best->name);
        updateStocksAfterProcess(best, stocks);
    }

    if (seq.size() < static_cast<size_t>(minSequenceLength) && !seq.empty()) {
        std::mt19937 &rng = randomGenerator;
        while (seq.size() < static_cast<size_t>(minSequenceLength)) seq.push_back(seq[rng() % seq.size()]);
    }

    return Individual(seq);
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

    int smartCount = populationSize * 0.8;

    for (int i = 0; i < smartCount; ++i) {
        population.push_back(createSmartIndividual());
    }

    for (int i = smartCount; i < populationSize; ++i) {
        population.push_back(createRandomIndividual());
    }

    cout << "Population initialized with " << populationSize
         << " individuals: " << smartCount << " smart, "
         << (populationSize - smartCount) << " random." << endl;
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
        selectedIndices.reserve(populationSize);
        uniform_int_distribution<size_t> randDist(0, population.size() - 1);
        for (int i = 0; i < populationSize; ++i) {
            selectedIndices.push_back(randDist(randomGenerator));
        }
        return selectedIndices;
    }

    if (minFitness == maxFitness) {
        selectedIndices.reserve(populationSize);
        uniform_int_distribution<size_t> randDist(0, population.size() - 1);
        for (int i = 0; i < populationSize; ++i) {
            selectedIndices.push_back(randDist(randomGenerator));
        }
        return selectedIndices;
    }

    double range = maxFitness - minFitness;
    vector<double> normalizedFitness;
    normalizedFitness.reserve(population.size());
    double totalNormalizedFitness = 0.0;

    for (const auto &individual : population) {
        double normalized = 0.001;
        if (isfinite(individual.fitnessScore)) {
            normalized = (individual.fitnessScore - minFitness) / range + 0.001;
        }
        normalizedFitness.push_back(normalized);
        totalNormalizedFitness += normalized;
    }

    if (totalNormalizedFitness <= 0.0 || !isfinite(totalNormalizedFitness)) {
        selectedIndices.reserve(populationSize);
        uniform_int_distribution<size_t> randDist(0, population.size() - 1);
        for (int i = 0; i < populationSize; ++i) {
            selectedIndices.push_back(randDist(randomGenerator));
        }
        return selectedIndices;
    }

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

    for (size_t i = 0; i < mutated.processSequence.size(); ++i) {
        if (probDist(randomGenerator) < mutationRate) {
            mutated.processSequence[i] =
                processNames[processDist(randomGenerator)];
        }
    }

    double structuralMutationRate = 0.02;
    if (probDist(randomGenerator) < structuralMutationRate &&
        mutated.processSequence.size() > 1) {
        uniform_int_distribution<size_t> posDist(
            0, mutated.processSequence.size());
        size_t pos = posDist(randomGenerator);

        if (probDist(randomGenerator) < 0.5 &&
            mutated.processSequence.size() >
                static_cast<size_t>(minSequenceLength)) {
            if (pos < mutated.processSequence.size()) {
                mutated.processSequence.erase(mutated.processSequence.begin() +
                                              pos);
            }
        } else {
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

        if (gen == 0 || currentBest.fitnessScore > bestOverall.fitnessScore) {
            bestOverall = currentBest;
            cout << "Generation " << (gen + 1) << "/" << generations
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