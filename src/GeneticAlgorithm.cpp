#include "GeneticAlgorithm.hpp"
#include <algorithm>
#include <iostream>
#include <numeric>
#include <unordered_map>
using namespace std;

namespace {
    double normalizeFitness(double score, double minF, double maxF) {
        return (isfinite(score) && maxF > minF)
                   ? (score - minF) / (maxF - minF) + 0.001
                   : 0.001;
    }
} // namespace

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

bool GeneticAlgorithm::canExecuteProcess(const Process *process,
                                         const map<string, int> &stocks) const {
    if (!process) return false;
    for (const auto &[resource, quantity] : process->inputs) {
        auto stockIt = stocks.find(resource);
        if (stockIt == stocks.end() || stockIt->second < quantity) return false;
    }
    return true;
}

void GeneticAlgorithm::updateStocksAfterProcess(
    const Process *process, map<string, int> &stocks) const {
    if (!process) return;
    for (const auto &[resource, quantity] : process->inputs)
        stocks[resource] -= quantity;
    for (const auto &[resource, quantity] : process->outputs)
        stocks[resource] += quantity;
}

Individual GeneticAlgorithm::createSmartIndividual() {
    map<string, int> stocks;
    for (const auto &stock : config.getStocks())
        stocks[stock.name] = stock.quantity;
    unordered_map<string, const Process *> processMap;
    for (const auto &process : config.getProcesses())
        processMap[process.name] = &process;
    vector<string> processSequence;
    int attempts = 0;
    while (processSequence.size() < (size_t)maxSequenceLength &&
           attempts < maxSequenceLength * 2) {
        ++attempts;
        vector<const Process *> executable;
        for (const auto &[name, ptr] : processMap)
            if (canExecuteProcess(ptr, stocks)) executable.push_back(ptr);
        if (executable.empty()) break;
        auto best = *min_element(
            executable.begin(), executable.end(),
            [&](const Process *a, const Process *b) {
                const auto &prio = simulator.getProcessPriority();
                int pa = prio.count(a->name) ? prio.at(a->name) : 3;
                int pb = prio.count(b->name) ? prio.at(b->name) : 3;
                return pa != pb ? pa < pb : a->nbCycle < b->nbCycle;
            });
        processSequence.push_back(best->name);
        updateStocksAfterProcess(best, stocks);
    }
    if (processSequence.size() < (size_t)minSequenceLength &&
        !processSequence.empty()) {
        std::mt19937 &rng = randomGenerator;
        while (processSequence.size() < (size_t)minSequenceLength)
            processSequence.push_back(
                processSequence[rng() % processSequence.size()]);
    }
    return Individual(processSequence);
}

Individual GeneticAlgorithm::createRandomIndividual() {
    uniform_int_distribution<> lengthDist(minSequenceLength, maxSequenceLength);
    uniform_int_distribution<> processDist(0, (int)processNames.size() - 1);
    int len = lengthDist(randomGenerator);
    vector<string> seq(len);
    generate(seq.begin(), seq.end(),
             [&]() { return processNames[processDist(randomGenerator)]; });
    return Individual(seq);
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
    for (auto &individual : population)
        individual.fitnessScore =
            simulator.runSimulation(individual.processSequence).fitness;
}

vector<size_t> GeneticAlgorithm::selectParents() {
    vector<size_t> selected;
    if (population.empty()) return selected;
    double minF = numeric_limits<double>::max(),
           maxF = numeric_limits<double>::lowest();
    for (const auto &ind : population)
        if (isfinite(ind.fitnessScore)) {
            minF = min(minF, ind.fitnessScore);
            maxF = max(maxF, ind.fitnessScore);
        }
    if (minF == numeric_limits<double>::max() ||
        maxF == numeric_limits<double>::lowest() || minF == maxF) {
        uniform_int_distribution<size_t> randDist(0, population.size() - 1);
        generate_n(back_inserter(selected), populationSize,
                   [&]() { return randDist(randomGenerator); });
        return selected;
    }
    vector<double> normF;
    double total = 0.0;
    for (const auto &ind : population) {
        double n = normalizeFitness(ind.fitnessScore, minF, maxF);
        normF.push_back(n);
        total += n;
    }
    uniform_real_distribution<> dist(0.0, total);
    for (int i = 0; i < populationSize; ++i) {
        double r = dist(randomGenerator), sum = 0.0;
        for (size_t j = 0; j < population.size(); ++j) {
            sum += normF[j];
            if (r <= sum) {
                selected.push_back(j);
                break;
            }
        }
        if (selected.size() <= i) selected.push_back(population.size() - 1);
    }
    return selected;
}

pair<Individual, Individual> GeneticAlgorithm::crossover(const Individual &p1,
                                                         const Individual &p2) {
    uniform_real_distribution<> crossDist(0.0, 1.0);
    if (crossDist(randomGenerator) > crossoverRate) return {p1, p2};
    const auto &s1 = p1.processSequence, &s2 = p2.processSequence;
    size_t l1 = s1.size(), l2 = s2.size(), shortLen = min(l1, l2);
    if (l1 < 2 || l2 < 2) return {p1, p2};
    uniform_int_distribution<size_t> pointDist(0, shortLen - 1);
    size_t pA = pointDist(randomGenerator), pB = pointDist(randomGenerator);
    if (pA == pB) pB = (pA + 1) % shortLen;
    if (pA > pB) swap(pA, pB);
    vector<string> c1, c2;
    c1.insert(c1.end(), s1.begin(), s1.begin() + pA);
    if (pB < l2)
        c1.insert(c1.end(), s2.begin() + pA, s2.begin() + pB);
    else if (pA < l2)
        c1.insert(c1.end(), s2.begin() + pA, s2.end());
    if (pB < l1) c1.insert(c1.end(), s1.begin() + pB, s1.end());
    c2.insert(c2.end(), s2.begin(), s2.begin() + pA);
    if (pB < l1)
        c2.insert(c2.end(), s1.begin() + pA, s1.begin() + pB);
    else if (pA < l1)
        c2.insert(c2.end(), s1.begin() + pA, s1.end());
    if (pB < l2) c2.insert(c2.end(), s2.begin() + pB, s2.end());
    return {Individual(c1), Individual(c2)};
}

Individual GeneticAlgorithm::mutate(const Individual &ind) {
    Individual m = ind;
    uniform_real_distribution<> probDist(0.0, 1.0);
    uniform_int_distribution<> processDist(0, (int)processNames.size() - 1);
    for (auto &name : m.processSequence)
        if (probDist(randomGenerator) < mutationRate)
            name = processNames[processDist(randomGenerator)];
    double structRate = 0.02;
    if (probDist(randomGenerator) < structRate &&
        m.processSequence.size() > 1) {
        uniform_int_distribution<size_t> posDist(0, m.processSequence.size());
        size_t pos = posDist(randomGenerator);
        if (probDist(randomGenerator) < 0.5 &&
            m.processSequence.size() > (size_t)minSequenceLength) {
            if (pos < m.processSequence.size())
                m.processSequence.erase(m.processSequence.begin() + pos);
        } else {
            string randomProcess = processNames[processDist(randomGenerator)];
            if (pos > m.processSequence.size()) pos = m.processSequence.size();
            m.processSequence.insert(m.processSequence.begin() + pos,
                                     randomProcess);
        }
    }
    return m;
}

void GeneticAlgorithm::selectNextGeneration() {
    vector<Individual> newPop;
    newPop.reserve(populationSize);
    sort(population.begin(), population.end());
    for (int i = 0; i < eliteCount && i < (int)population.size(); ++i)
        newPop.push_back(population[i]);
    vector<size_t> parentIdx = selectParents();
    if (parentIdx.empty() && populationSize > eliteCount) {
        while (newPop.size() < (size_t)populationSize)
            newPop.push_back(!population.empty() ? population[0]
                                                 : createSmartIndividual());
        population = newPop;
        return;
    }
    if (parentIdx.empty()) {
        population = newPop;
        return;
    }
    uniform_int_distribution<size_t> parentDist(0, parentIdx.size() - 1);
    while (newPop.size() < (size_t)populationSize) {
        size_t idx1 = parentDist(randomGenerator),
               idx2 = parentDist(randomGenerator);
        if (idx1 >= parentIdx.size() || idx2 >= parentIdx.size()) continue;
        size_t p1 = parentIdx[idx1], p2 = parentIdx[idx2];
        if (p1 >= population.size() || p2 >= population.size()) continue;
        int tries = 0;
        while (p1 == p2 && parentIdx.size() > 1 && tries < 10) {
            idx2 = parentDist(randomGenerator);
            if (idx2 >= parentIdx.size()) continue;
            p2 = parentIdx[idx2];
            if (p2 >= population.size()) continue;
            tries++;
        }
        auto [c1, c2] = crossover(population[p1], population[p2]);
        c1 = mutate(c1);
        c2 = mutate(c2);
        if (newPop.size() < (size_t)populationSize) newPop.push_back(c1);
        if (newPop.size() < (size_t)populationSize) newPop.push_back(c2);
    }
    population = newPop;
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
    if (finalBest.fitnessScore > bestOverall.fitnessScore)
        bestOverall = finalBest;
    cout << "Evolution finished. Final best fitness: "
         << bestOverall.fitnessScore << endl;
    return bestOverall;
}