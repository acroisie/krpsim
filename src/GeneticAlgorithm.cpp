#include "GeneticAlgorithm.hpp"
#include <algorithm>
#include <future>
#include <iostream>
#include <numeric>
#include <thread>
#include <unordered_map>
using namespace std;

// Utility type for clarity
using Stocks = std::map<std::string, int>;

namespace {
    // Random selection of indices in the population
    vector<size_t> randomSelection(size_t count, size_t populationSize,
                                   mt19937 &rng) {
        vector<size_t> indices;
        uniform_int_distribution<size_t> randDist(0, populationSize - 1);
        for (size_t i = 0; i < count; ++i) indices.push_back(randDist(rng));
        return indices;
    }
} // namespace

GeneticAlgorithm::GeneticAlgorithm(const Config &cfg, Simulator &sim,
                                   int popSize, double mutRate,
                                   double crossRate, int elite, int minLen,
                                   int maxLen)
    : config(cfg), simulator(sim), populationSize(popSize),
      mutationRate(mutRate), crossoverRate(crossRate),
      eliteCount(min(max(0, elite), popSize)), minSequenceLength(minLen),
      maxSequenceLength(maxLen), randomGenerator(random_device{}()) {
    for (const auto &p : config.getProcesses()) processNames.push_back(p.name);
}

// Check if a process can be started with the available stocks
bool GeneticAlgorithm::canStartProcess(const Process *process,
                                       const Stocks &availableStocks) const {
    if (!process) return false;
    for (const auto &[resource, quantity] : process->inputs)
        if (availableStocks.count(resource) == 0 ||
            availableStocks.at(resource) < quantity)
            return false;
    return true;
}

// Update stocks after a process execution
void GeneticAlgorithm::updateStocksAfterProcess(const Process *process,
                                                Stocks &stocks) const {
    if (!process) return;
    for (const auto &[resource, quantity] : process->inputs)
        stocks[resource] -= quantity;
    for (const auto &[resource, quantity] : process->outputs)
        stocks[resource] += quantity;
}

// Create a smart individual by always picking the best process available
Individual GeneticAlgorithm::createSmartIndividual() {
    Stocks stocks;
    for (const auto &stock : config.getStocks())
        stocks[stock.name] = stock.quantity;
    std::unordered_map<std::string, const Process *> processByName;
    for (const auto &process : config.getProcesses())
        processByName[process.name] = &process;
    std::vector<std::string> sequence;
    int attempts = 0;
    while (sequence.size() < static_cast<size_t>(maxSequenceLength) &&
           attempts++ < maxSequenceLength * 2) {
        std::vector<const Process *> availableProcesses;
        for (const auto &[name, ptr] : processByName)
            if (canStartProcess(ptr, stocks)) availableProcesses.push_back(ptr);
        if (availableProcesses.empty()) break;
        // Pick the process with the highest priority and shortest duration
        auto best = *std::min_element(
            availableProcesses.begin(), availableProcesses.end(),
            [&](const Process *a, const Process *b) {
                const auto &priority = simulator.getProcessPriority();
                int pa = priority.count(a->name) ? priority.at(a->name) : 3;
                int pb = priority.count(b->name) ? priority.at(b->name) : 3;
                return pa != pb ? pa < pb : a->nbCycle < b->nbCycle;
            });
        sequence.push_back(best->name);
        updateStocksAfterProcess(best, stocks);
    }
    // If the sequence is too short, repeat useful choices
    if (sequence.size() < static_cast<size_t>(minSequenceLength) &&
        !sequence.empty()) {
        std::mt19937 &rng = randomGenerator;
        while (sequence.size() < static_cast<size_t>(minSequenceLength))
            sequence.push_back(sequence[rng() % sequence.size()]);
    }
    return Individual(sequence);
}

Individual GeneticAlgorithm::createRandomIndividual() {
    uniform_int_distribution<> lenDist(minSequenceLength, maxSequenceLength);
    uniform_int_distribution<> procDist(0, (int)processNames.size() - 1);
    int len = lenDist(randomGenerator);
    vector<string> sequence(len);
    for (int i = 0; i < len; ++i)
        sequence[i] = processNames[procDist(randomGenerator)];
    return Individual(sequence);
}

void GeneticAlgorithm::initializePopulation() {
    population.clear();
    population.reserve(populationSize);
    int smartCount = populationSize * 0.8;
    for (int i = 0; i < smartCount; ++i)
        population.push_back(createSmartIndividual());
    for (int i = smartCount; i < populationSize; ++i)
        population.push_back(createRandomIndividual());
    cout << "Population initialized with " << populationSize
         << " individuals: " << smartCount << " smart, "
         << (populationSize - smartCount) << " random." << endl;
}

void GeneticAlgorithm::evaluatePopulation() {
    // Parallel evaluation of the population
    size_t numThreads = std::thread::hardware_concurrency();
    if (numThreads == 0) numThreads = 4; // fallback
    size_t popSize = population.size();
    std::vector<std::future<void>> futures;
    size_t chunk = (popSize + numThreads - 1) / numThreads;
    for (size_t t = 0; t < numThreads; ++t) {
        size_t start = t * chunk;
        size_t end = std::min(start + chunk, popSize);
        if (start >= end) break;
        futures.push_back(std::async(std::launch::async, [this, start, end]() {
            for (size_t i = start; i < end; ++i) {
                population[i].fitnessScore =
                    simulator.runSimulation(population[i].processSequence)
                        .fitness;
            }
        }));
    }
    for (auto &f : futures) f.get();
}

vector<size_t> GeneticAlgorithm::selectParents() {
    if (population.empty()) return {};
    double minFitness = numeric_limits<double>::max(),
           maxFitness = numeric_limits<double>::lowest();
    for (const auto &individual : population)
        if (isfinite(individual.fitnessScore)) {
            minFitness = min(minFitness, individual.fitnessScore);
            maxFitness = max(maxFitness, individual.fitnessScore);
        }
    if (minFitness == numeric_limits<double>::max() ||
        maxFitness == numeric_limits<double>::lowest() ||
        minFitness == maxFitness) {
        return randomSelection(populationSize, population.size(),
                               randomGenerator);
    }
    double range = maxFitness - minFitness;
    vector<double> normalizedFitness;
    normalizedFitness.reserve(population.size());
    double total = 0.0;
    for (const auto &individual : population) {
        double normalized = 0.001;
        if (isfinite(individual.fitnessScore))
            normalized = (individual.fitnessScore - minFitness) / range + 0.001;
        normalizedFitness.push_back(normalized);
        total += normalized;
    }
    vector<size_t> selectedIndices;
    uniform_real_distribution<> dist(0.0, total);
    for (int i = 0; i < populationSize; ++i) {
        double randomPoint = dist(randomGenerator), sum = 0.0;
        for (size_t j = 0; j < population.size(); ++j) {
            sum += normalizedFitness[j];
            if (randomPoint <= sum) {
                selectedIndices.push_back(j);
                break;
            }
        }
        if (selectedIndices.size() < static_cast<size_t>(i + 1) &&
            !population.empty())
            selectedIndices.push_back(population.size() - 1);
    }
    return selectedIndices;
}

pair<Individual, Individual>
GeneticAlgorithm::crossover(const Individual &parent1,
                            const Individual &parent2) {
    uniform_real_distribution<> crossDist(0.0, 1.0);
    if (crossDist(randomGenerator) > crossoverRate) return {parent1, parent2};
    const auto &seq1 = parent1.processSequence, &seq2 = parent2.processSequence;
    size_t len1 = seq1.size(), len2 = seq2.size(), shortLen = min(len1, len2);
    if (len1 < 2 || len2 < 2) return {parent1, parent2};
    uniform_int_distribution<size_t> pointDist(0, shortLen - 1);
    size_t pointA = pointDist(randomGenerator),
           pointB = pointDist(randomGenerator);
    if (pointA == pointB) pointB = (pointA + 1) % shortLen;
    if (pointA > pointB) swap(pointA, pointB);
    vector<string> childSeq1, childSeq2;
    childSeq1.insert(childSeq1.end(), seq1.begin(), seq1.begin() + pointA);
    if (pointB < len2)
        childSeq1.insert(childSeq1.end(), seq2.begin() + pointA,
                         seq2.begin() + pointB);
    else if (pointA < len2)
        childSeq1.insert(childSeq1.end(), seq2.begin() + pointA, seq2.end());
    if (pointB < len1)
        childSeq1.insert(childSeq1.end(), seq1.begin() + pointB, seq1.end());
    childSeq2.insert(childSeq2.end(), seq2.begin(), seq2.begin() + pointA);
    if (pointB < len1)
        childSeq2.insert(childSeq2.end(), seq1.begin() + pointA,
                         seq1.begin() + pointB);
    else if (pointA < len1)
        childSeq2.insert(childSeq2.end(), seq1.begin() + pointA, seq1.end());
    if (pointB < len2)
        childSeq2.insert(childSeq2.end(), seq2.begin() + pointB, seq2.end());
    return {Individual(childSeq1), Individual(childSeq2)};
}

Individual GeneticAlgorithm::mutate(const Individual &individual) {
    Individual mutated = individual;
    uniform_real_distribution<> probDist(0.0, 1.0);
    uniform_int_distribution<> procDist(0, (int)processNames.size() - 1);
    for (size_t i = 0; i < mutated.processSequence.size(); ++i)
        if (probDist(randomGenerator) < mutationRate)
            mutated.processSequence[i] =
                processNames[procDist(randomGenerator)];
    double structRate = 0.02;
    if (probDist(randomGenerator) < structRate &&
        mutated.processSequence.size() > 1) {
        uniform_int_distribution<size_t> posDist(
            0, mutated.processSequence.size());
        size_t pos = posDist(randomGenerator);
        if (probDist(randomGenerator) < 0.5 &&
            mutated.processSequence.size() > (size_t)minSequenceLength) {
            if (pos < mutated.processSequence.size())
                mutated.processSequence.erase(mutated.processSequence.begin() +
                                              pos);
        } else {
            string randomProcess = processNames[procDist(randomGenerator)];
            if (pos > mutated.processSequence.size())
                pos = mutated.processSequence.size();
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
    for (int i = 0; i < eliteCount && i < (int)population.size(); ++i)
        newPopulation.push_back(population[i]);
    vector<size_t> parentIndices = selectParents();
    if (parentIndices.empty() && populationSize > eliteCount) {
        while (newPopulation.size() < (size_t)populationSize)
            newPopulation.push_back(population.empty() ? createSmartIndividual()
                                                       : population[0]);
        population = newPopulation;
        return;
    }
    if (parentIndices.empty()) {
        population = newPopulation;
        return;
    }
    uniform_int_distribution<size_t> parentDist(0, parentIndices.size() - 1);
    while (newPopulation.size() < (size_t)populationSize) {
        size_t idx1 = parentDist(randomGenerator),
               idx2 = parentDist(randomGenerator);
        if (idx1 >= parentIndices.size() || idx2 >= parentIndices.size())
            continue;
        size_t parent1Idx = parentIndices[idx1],
               parent2Idx = parentIndices[idx2];
        if (parent1Idx >= population.size() || parent2Idx >= population.size())
            continue;
        int tries = 0;
        while (parent1Idx == parent2Idx && parentIndices.size() > 1 &&
               tries++ < 10) {
            idx2 = parentDist(randomGenerator);
            if (idx2 >= parentIndices.size()) continue;
            parent2Idx = parentIndices[idx2];
            if (parent2Idx >= population.size()) continue;
        }
        auto [child1, child2] =
            crossover(population[parent1Idx], population[parent2Idx]);
        child1 = mutate(child1);
        child2 = mutate(child2);
        if (newPopulation.size() < (size_t)populationSize)
            newPopulation.push_back(child1);
        if (newPopulation.size() < (size_t)populationSize)
            newPopulation.push_back(child2);
    }
    population = newPopulation;
}

Individual GeneticAlgorithm::getBestIndividual() const {
    if (population.empty()) return Individual();
    return *max_element(population.begin(), population.end(),
                        [](const Individual &a, const Individual &b) {
                            return a.fitnessScore < b.fitnessScore;
                        });
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
    if (finalBest.fitnessScore > bestOverall.fitnessScore)
        bestOverall = finalBest;
    cout << "Evolution finished. Final best fitness: "
         << bestOverall.fitnessScore << endl;
    return bestOverall;
}