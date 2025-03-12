#include "ProcessManager.hpp"
#include "Config.hpp"
#include "Individual.hpp"
#include "Process.hpp"
#include <algorithm>
#include <iostream>
#include <limits>
#include <random>

using namespace std;

const int ProcessManager::POPULATION_SIZE = 20;

ProcessManager::ProcessManager(const Config &config, int delayLimit)
    : config_(config), delayLimit_(delayLimit), currentCycle_(0) {
    for (const auto &stock : config_.getStocks()) {
        currentStocks_[stock.name] = stock.quantity;
    }
}

void ProcessManager::initializePopulation() {
    population_.resize(POPULATION_SIZE);
    vector<string> allProcessNames;
    for (const auto &proc : config_.getProcesses()) {
        allProcessNames.push_back(proc.name);
    }

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> lengthDist(1, delayLimit_);
    uniform_int_distribution<> procDist(0, (int)allProcessNames.size() - 1);

    for (int i = 0; i < POPULATION_SIZE; ++i) {
        int seqLen = lengthDist(gen);
        vector<string> sequence;
        sequence.reserve(seqLen);

        for (int j = 0; j < seqLen; ++j) {
            sequence.push_back(allProcessNames[procDist(gen)]);
        }
        population_[i] = Individual(sequence);
    }
}

void ProcessManager::evaluateFitness() {
    for (auto &indiv : population_) {
        indiv.fitness = calculateFitness(indiv);
    }
}

double ProcessManager::calculateFitness(Individual &individual) {
    currentStocks_.clear();
    for (const auto &stock : config_.getStocks()) {
        currentStocks_[stock.name] = stock.quantity;
    }
    executionLogs_.clear();
    currentCycle_ = 0;

    struct RunningProcess {
        const Process *process;
        int completionCycle;
    };
    vector<RunningProcess> runningProcesses;

    size_t index = 0;
    while (currentCycle_ < delayLimit_ && index < individual.processSequence.size()) {
        // Libérer les process terminés
        auto it = runningProcesses.begin();
        while (it != runningProcesses.end()) {
            if (it->completionCycle <= currentCycle_) {
                // On applique les outputs
                for (const auto &out : it->process->outputs) {
                    currentStocks_[out.first] += out.second;
                }
                it = runningProcesses.erase(it);
            } else {
                ++it;
            }
        }

        // Tenter d’exécuter le process courant
        const string &procName = individual.processSequence[index];
        vector<const Process *> runnable = getRunnableProcesses();
        const Process *chosen = nullptr;

        for (auto p : runnable) {
            if (p->name == procName) {
                chosen = p;
                break;
            }
        }

        if (chosen) {
            // Consommer les inputs
            for (const auto &in : chosen->inputs) {
                currentStocks_[in.first] -= in.second;
            }
            runningProcesses.push_back({chosen, currentCycle_ + chosen->nbCycle});
            executionLogs_.push_back({currentCycle_, chosen->name});
            ++index;
        } else {
            // Si on ne peut pas le lancer
            if (!runningProcesses.empty()) {
                int nextCompletion = delayLimit_;
                for (auto &rp : runningProcesses) {
                    nextCompletion = min(nextCompletion, rp.completionCycle);
                }
                currentCycle_ = nextCompletion;
            } else {
                // Personne ne tourne, on avance
                ++index;
                ++currentCycle_;
            }
        }
    }

    // Laisser finir les process en cours
    while (!runningProcesses.empty() && currentCycle_ < delayLimit_) {
        int nextCompletion = delayLimit_;
        for (auto &rp : runningProcesses) {
            nextCompletion = min(nextCompletion, rp.completionCycle);
        }
        currentCycle_ = nextCompletion;

        auto it = runningProcesses.begin();
        while (it != runningProcesses.end()) {
            if (it->completionCycle <= currentCycle_) {
                for (const auto &out : it->process->outputs) {
                    currentStocks_[out.first] += out.second;
                }
                it = runningProcesses.erase(it);
            } else {
                ++it;
            }
        }
    }

    double fitness = 0.0;
    const auto &goals = config_.getOptimizeGoal();
    if (!goals.empty()) {
        // On gère un seul objectif principal ou le premier trouvé
        std::string goal = goals[0];
        if (goal == "time") {
            fitness = -static_cast<double>(currentCycle_);
        } else {
            // Maximiser la quantité d’un stock
            if (currentStocks_.count(goal)) {
                fitness = currentStocks_[goal];
            } else {
                fitness = numeric_limits<double>::lowest();
            }
        }
    }
    return fitness;
}

vector<const Process *> ProcessManager::getRunnableProcesses() {
    vector<const Process *> result;
    for (const auto &proc : config_.getProcesses()) {
        bool canRun = true;
        for (const auto &in : proc.inputs) {
            if (currentStocks_[in.first] < in.second) {
                canRun = false;
                break;
            }
        }
        if (canRun) {
            result.push_back(&proc);
        }
    }
    return result;
}

bool ProcessManager::executeProcess(const Process *process) {
    if (!process) return false;
    for (const auto &in : process->inputs) {
        currentStocks_[in.first] -= in.second;
        if (currentStocks_[in.first] < 0) {
            return false;
        }
    }
    for (const auto &out : process->outputs) {
        currentStocks_[out.first] += out.second;
    }
    return true;
}

vector<Individual> ProcessManager::selectParents(const vector<Individual> &population) {
    vector<Individual> parents;
    if (population.empty()) {
        return parents;
    }

    double totalFitness = 0.0;
    double minFitness = numeric_limits<double>::max();

    for (const auto &indiv : population) {
        totalFitness += indiv.fitness;
        if (indiv.fitness < minFitness) {
            minFitness = indiv.fitness;
        }
    }

    // On shift si jamais toutes les fitness sont négatives
    vector<double> adjustedFitness;
    adjustedFitness.reserve(population.size());

    if (totalFitness <= 0.0) {
        double shift = std::abs(minFitness) + 1.0;
        totalFitness = 0.0;
        for (const auto &ind : population) {
            double adj = ind.fitness + shift;
            adjustedFitness.push_back(adj);
            totalFitness += adj;
        }
    } else {
        for (const auto &ind : population) {
            adjustedFitness.push_back(ind.fitness);
        }
    }

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dist(0.0, totalFitness);

    int nbParents = min(POPULATION_SIZE, (int)population.size());
    for (int i = 0; i < nbParents; ++i) {
        double r = dist(gen);
        double accum = 0.0;
        for (size_t j = 0; j < population.size(); ++j) {
            accum += adjustedFitness[j];
            if (r <= accum) {
                parents.push_back(population[j]);
                break;
            }
        }
        if ((int)parents.size() <= i && !population.empty()) {
            parents.push_back(population[0]);
        }
    }
    return parents;
}

pair<Individual, Individual> ProcessManager::crossover(const Individual &parent1, const Individual &parent2) {
    if (parent1.processSequence.empty() || parent2.processSequence.empty()) {
        return make_pair(parent1, parent2);
    }
    random_device rd;
    mt19937 gen(rd());
    size_t minSize = min(parent1.processSequence.size(), parent2.processSequence.size());
    if (minSize == 0) {
        return make_pair(parent1, parent2);
    }

    uniform_int_distribution<> dist(0, (int)minSize - 1);
    size_t crossoverPoint = dist(gen);

    vector<string> child1, child2;
    child1.reserve(parent1.processSequence.size() + parent2.processSequence.size());
    child2.reserve(parent1.processSequence.size() + parent2.processSequence.size());

    for (size_t i = 0; i < crossoverPoint; ++i) {
        if (i < parent1.processSequence.size()) {
            child1.push_back(parent1.processSequence[i]);
        }
        if (i < parent2.processSequence.size()) {
            child2.push_back(parent2.processSequence[i]);
        }
    }
    for (size_t i = crossoverPoint; i < parent2.processSequence.size(); ++i) {
        child1.push_back(parent2.processSequence[i]);
    }
    for (size_t i = crossoverPoint; i < parent1.processSequence.size(); ++i) {
        child2.push_back(parent1.processSequence[i]);
    }

    return make_pair(Individual(child1), Individual(child2));
}

Individual ProcessManager::mutate(const Individual &individual, double mutationRate) {
    auto allProcs = config_.getProcesses();
    if (allProcs.empty() || individual.processSequence.empty()) {
        return individual;
    }

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> distProb(0.0, 1.0);
    uniform_int_distribution<> distProc(0, (int)allProcs.size() - 1);

    vector<string> mutatedSeq = individual.processSequence;
    for (size_t i = 0; i < mutatedSeq.size(); ++i) {
        if (distProb(gen) < mutationRate) {
            mutatedSeq[i] = allProcs[distProc(gen)].name;
        }
    }
    return Individual(mutatedSeq);
}

void ProcessManager::runGeneticAlgorithm() {
    initializePopulation();
    const int NUM_GENERATIONS = 20;
    const double MUTATION_RATE = 0.1;

    evaluateFitness();

    for (int generation = 0; generation < NUM_GENERATIONS; ++generation) {
        vector<Individual> parents = selectParents(population_);
        vector<Individual> newPopulation;

        if (!population_.empty()) {
            auto bestInPop = findBestIndividual(population_);
            newPopulation.push_back(bestInPop);
        }

        while (newPopulation.size() < population_.size() && parents.size() >= 2) {
            size_t idx1 = rand() % parents.size();
            size_t idx2 = rand() % parents.size();
            while (idx2 == idx1 && parents.size() > 1) {
                idx2 = rand() % parents.size();
            }

            auto children = crossover(parents[idx1], parents[idx2]);
            auto child1 = mutate(children.first, MUTATION_RATE);
            auto child2 = mutate(children.second, MUTATION_RATE);

            if (newPopulation.size() < population_.size()) {
                newPopulation.push_back(child1);
            }
            if (newPopulation.size() < population_.size()) {
                newPopulation.push_back(child2);
            }
        }

        if (!newPopulation.empty()) {
            population_ = newPopulation;
        }
        evaluateFitness();
    }

    if (!population_.empty()) {
        bestSolution_ = findBestIndividual(population_);
        generateOutput();
    } else {
        cout << "Error: Population is empty after evolution!" << endl;
    }
}

Individual ProcessManager::findBestIndividual(const vector<Individual> &population) {
    Individual best = population.front();
    for (const auto &ind : population) {
        if (ind.fitness > best.fitness) {
            best = ind;
        }
    }
    return best;
}

void ProcessManager::generateOutput() {
    // On reprend la simulation avec bestSolution_ pour afficher les logs
    // On reset les stocks
    currentStocks_.clear();
    for (const auto &stock : config_.getStocks()) {
        currentStocks_[stock.name] = stock.quantity;
    }
    executionLogs_.clear();
    currentCycle_ = 0;

    struct RunningProcess {
        const Process *process;
        int completionCycle;
    };
    vector<RunningProcess> runningProcesses;

    const auto &sequence = bestSolution_.processSequence;
    size_t idx = 0;

    while (currentCycle_ < delayLimit_ && idx < sequence.size()) {
        auto it = runningProcesses.begin();
        while (it != runningProcesses.end()) {
            if (it->completionCycle <= currentCycle_) {
                for (auto &out : it->process->outputs) {
                    currentStocks_[out.first] += out.second;
                }
                it = runningProcesses.erase(it);
            } else {
                ++it;
            }
        }

        const string &procName = sequence[idx];
        vector<const Process *> runnable = getRunnableProcesses();
        const Process *chosen = nullptr;

        for (auto &p : runnable) {
            if (p->name == procName) {
                chosen = p;
                break;
            }
        }

        if (chosen) {
            for (auto &in : chosen->inputs) {
                currentStocks_[in.first] -= in.second;
            }
            runningProcesses.push_back({chosen, currentCycle_ + chosen->nbCycle});
            executionLogs_.push_back({currentCycle_, chosen->name});
            ++idx;
        } else {
            if (!runningProcesses.empty()) {
                int nextCompletion = delayLimit_;
                for (auto &rp : runningProcesses) {
                    nextCompletion = min(nextCompletion, rp.completionCycle);
                }
                currentCycle_ = nextCompletion;
            } else {
                ++idx;
                ++currentCycle_;
            }
        }
    }

    while (!runningProcesses.empty() && currentCycle_ < delayLimit_) {
        int nextCompletion = delayLimit_;
        for (auto &rp : runningProcesses) {
            nextCompletion = min(nextCompletion, rp.completionCycle);
        }
        currentCycle_ = nextCompletion;

        auto it = runningProcesses.begin();
        while (it != runningProcesses.end()) {
            if (it->completionCycle <= currentCycle_) {
                for (auto &out : it->process->outputs) {
                    currentStocks_[out.first] += out.second;
                }
                it = runningProcesses.erase(it);
            } else {
                ++it;
            }
        }
    }

    cout << "Nice file! " << config_.getProcesses().size() << " processes, " << config_.getStocks().size()
         << " stocks, " << config_.getOptimizeGoal().size() << " to optimize" << endl;

    cout << "Main walk:" << endl;
    for (auto &log : executionLogs_) {
        cout << log.first << ":" << log.second << endl;
    }

    cout << "no more process doable at time " << currentCycle_ << endl;

    cout << "Stock :" << endl;
    for (auto &stock : currentStocks_) {
        cout << stock.first << "=> " << stock.second << endl;
    }
}
