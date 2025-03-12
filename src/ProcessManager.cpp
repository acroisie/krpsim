#include "ProcessManager.hpp"
#include "Config.hpp"
#include "Individual.hpp"
#include "Process.hpp"
#include <algorithm>
#include <chrono>
#include <chrono> // Ajouter cet include en haut du fichier
#include <iostream>
#include <limits>
#include <random>
#include <set>

using namespace std;

const int ProcessManager::POPULATION_SIZE = 500;   // Population plus large
const int ProcessManager::MAX_GENERATIONS = 300;   // Plus de générations
const double ProcessManager::MUTATION_RATE = 0.25; // Mutation plus agressive
const double ProcessManager::ELITISM_RATE = 0.15;  // Plus d'élitisme
const double ProcessManager::CROSSOVER_RATE = 0.9; // Croisement plus fréquent
const int ProcessManager::TOURNAMENT_SIZE = 6;     // Tournoi plus sélectif

ProcessManager::ProcessManager(const Config &config, int delayLimit)
    : config_(config), delayLimit_(delayLimit), currentCycle_(0) {
    for (const auto &stock : config_.getStocks()) {
        currentStocks_[stock.name] = stock.quantity;
    }
}

void ProcessManager::initializePopulation() {
    population_.resize(POPULATION_SIZE);
    vector<string> allProcessNames;
    vector<string> buyProcesses;  // Processus d'achat
    vector<string> prodProcesses; // Processus de production
    vector<string> sellProcesses; // Processus de vente

    for (const auto &proc : config_.getProcesses()) {
        allProcessNames.push_back(proc.name);

        // Catégoriser les processus
        if (proc.name.find("buy_") != string::npos) {
            buyProcesses.push_back(proc.name);
        } else if (proc.name.find("vente_") != string::npos) {
            sellProcesses.push_back(proc.name);
        } else if (proc.name.find("do_") != string::npos) {
            prodProcesses.push_back(proc.name);
        }
    }

    random_device rd;
    mt19937 gen(rd());

    // Séquences plus longues
    uniform_int_distribution<> lengthDist(delayLimit_ / 5, delayLimit_ / 2);

    // Répartir les stratégies sur toute la population
    for (int i = 0; i < POPULATION_SIZE; ++i) {
        vector<string> sequence;

        if (i < POPULATION_SIZE * 0.3 && !buyProcesses.empty() && !prodProcesses.empty() && !sellProcesses.empty()) {
            // Stratégie intelligente: acheter -> produire -> vendre
            int cyclesRequired = delayLimit_ / 3;

            // Phase 1: Acheter beaucoup de matières premières
            for (int j = 0; j < cyclesRequired / 3; ++j) {
                for (const auto &buy : buyProcesses) {
                    sequence.push_back(buy);
                }
            }

            // Phase 2: Production
            for (int j = 0; j < cyclesRequired / 3; ++j) {
                for (const auto &prod : prodProcesses) {
                    sequence.push_back(prod);
                }
            }

            // Phase 3: Vendre
            for (int j = 0; j < cyclesRequired / 3; ++j) {
                for (const auto &sell : sellProcesses) {
                    sequence.push_back(sell);
                }
            }

            // Répéter le cycle acheter -> produire -> vendre
            int remainingLength = lengthDist(gen) - sequence.size();
            for (int j = 0; j < remainingLength; ++j) {
                if (j % 3 == 0 && !buyProcesses.empty()) {
                    sequence.push_back(buyProcesses[gen() % buyProcesses.size()]);
                } else if (j % 3 == 1 && !prodProcesses.empty()) {
                    sequence.push_back(prodProcesses[gen() % prodProcesses.size()]);
                } else if (!sellProcesses.empty()) {
                    sequence.push_back(sellProcesses[gen() % sellProcesses.size()]);
                }
            }
        } else if (i < POPULATION_SIZE * 0.5) {
            // Séquences variées avec beaucoup de répétitions
            int seqLen = lengthDist(gen);
            for (size_t j = 0; j < 10 && j < allProcessNames.size(); ++j) {
                string proc = allProcessNames[gen() % allProcessNames.size()];
                for (int k = 0; k < seqLen / 10; ++k) {
                    sequence.push_back(proc);
                }
            }
        } else if (i < POPULATION_SIZE * 0.7) {
            // Séquences avec alternance
            int patternSize = 1 + gen() % 3; // Patterns de taille 1-3
            vector<string> pattern;
            for (int j = 0; j < patternSize; ++j) {
                pattern.push_back(allProcessNames[gen() % allProcessNames.size()]);
            }

            int seqLen = lengthDist(gen);
            for (int j = 0; j < seqLen; ++j) {
                sequence.push_back(pattern[j % patternSize]);
            }
        } else {
            // Séquence aléatoire pour le reste
            int seqLen = lengthDist(gen);
            sequence.reserve(seqLen);
            for (int j = 0; j < seqLen; ++j) {
                sequence.push_back(allProcessNames[gen() % allProcessNames.size()]);
            }
        }

        population_[i] = Individual(sequence);
    }
}

void ProcessManager::evaluateFitness() {
#pragma omp parallel for
    for (size_t i = 0; i < population_.size(); ++i) {
        population_[i].fitness = calculateFitness(population_[i]);
    }
}

double ProcessManager::calculateFitness(Individual &individual) {
    // Limiter le nombre maximal de processus à évaluer
    const size_t MAX_PROCESSES_TO_EVALUATE = 1000;
    if (individual.processSequence.size() > MAX_PROCESSES_TO_EVALUATE) {
        individual.processSequence.resize(MAX_PROCESSES_TO_EVALUATE);
    }

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
        std::string goal = goals[0];
        if (goal == "time") {
            fitness = -static_cast<double>(currentCycle_);
        } else {
            // Objectif principal (euros ou autre)
            if (currentStocks_.count(goal)) {
                fitness = currentStocks_[goal];

                // Bonus pour l'activité générale - applicable à n'importe quel problème
                double activityBonus = 0.0;

                // 1. Bonus pour nombre de processus différents exécutés
                std::set<std::string> uniqueProcesses;
                for (const auto &log : executionLogs_) {
                    uniqueProcesses.insert(log.second);
                }
                activityBonus += uniqueProcesses.size() * 5.0;

                // 2. Bonus pour utilisation des ressources
                size_t resourceTypes = 0;
                double resourceUtilization = 0.0;
                for (const auto &stock : currentStocks_) {
                    if (stock.second > 0 && stock.first != goal) {
                        resourceTypes++;
                        resourceUtilization += stock.second;
                    }
                }
                // Utiliser resourceUtilization dans le calcul du bonus
                activityBonus += resourceTypes * 10.0 + std::min(resourceUtilization * 0.01, 50.0);

                // 3. Bonus pour longueur d'exécution
                activityBonus += std::min(currentCycle_, 100) * 2.0;

                // Limiter le bonus à 10% de la fitness principale
                fitness += std::min(activityBonus, fitness * 0.1);
            } else {
                fitness = std::numeric_limits<double>::lowest();
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

    // Utiliser deux points de croisement plutôt qu'un seul
    vector<string> child1, child2;
    size_t size1 = parent1.processSequence.size();
    size_t size2 = parent2.processSequence.size();

    if (size1 < 2 || size2 < 2) {
        return make_pair(parent1, parent2);
    }

    uniform_int_distribution<> dist1(0, size1 - 1);
    uniform_int_distribution<> dist2(0, size2 - 1);

    size_t point1_p1 = dist1(gen);
    size_t point2_p1 = dist1(gen);
    if (point1_p1 > point2_p1) std::swap(point1_p1, point2_p1);

    size_t point1_p2 = dist2(gen);
    size_t point2_p2 = dist2(gen);
    if (point1_p2 > point2_p2) std::swap(point1_p2, point2_p2);

    // Construire enfant 1: début parent1 + milieu parent2 + fin parent1
    child1.insert(child1.end(), parent1.processSequence.begin(), parent1.processSequence.begin() + point1_p1);
    child1.insert(child1.end(), parent2.processSequence.begin() + point1_p2,
                  parent2.processSequence.begin() + point2_p2);
    child1.insert(child1.end(), parent1.processSequence.begin() + point2_p1, parent1.processSequence.end());

    // Construire enfant 2: début parent2 + milieu parent1 + fin parent2
    child2.insert(child2.end(), parent2.processSequence.begin(), parent2.processSequence.begin() + point1_p2);
    child2.insert(child2.end(), parent1.processSequence.begin() + point1_p1,
                  parent1.processSequence.begin() + point2_p1);
    child2.insert(child2.end(), parent2.processSequence.begin() + point2_p2, parent2.processSequence.end());

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

    // 20% de chance d'ajouter des processus à la fin (au lieu de 10%)
    if (distProb(gen) < 0.2 && mutatedSeq.size() < static_cast<size_t>(delayLimit_ / 5)) {
        int addCount = 5 + gen() % 20; // Ajouter entre 5 et 25 processus (plus qu'avant)

        // Ajouter des séquences intelligentes lors de la mutation
        vector<string> buyProcs, prodProcs, sellProcs;
        for (const auto &proc : allProcs) {
            if (proc.name.find("buy_") != string::npos) {
                buyProcs.push_back(proc.name);
            } else if (proc.name.find("vente_") != string::npos) {
                sellProcs.push_back(proc.name);
            } else if (proc.name.find("do_") != string::npos) {
                prodProcs.push_back(proc.name);
            }
        }

        // 50% de chance d'ajouter une séquence intelligente
        if (distProb(gen) < 0.5 && !buyProcs.empty() && !prodProcs.empty() && !sellProcs.empty()) {
            for (int i = 0; i < addCount; ++i) {
                if (i % 3 == 0) {
                    mutatedSeq.push_back(buyProcs[gen() % buyProcs.size()]);
                } else if (i % 3 == 1) {
                    mutatedSeq.push_back(prodProcs[gen() % prodProcs.size()]);
                } else {
                    mutatedSeq.push_back(sellProcs[gen() % sellProcs.size()]);
                }
            }
        } else {
            for (int i = 0; i < addCount; ++i) {
                mutatedSeq.push_back(allProcs[distProc(gen)].name);
            }
        }
    }

    // Mutation standard des processus existants
    for (size_t i = 0; i < mutatedSeq.size(); ++i) {
        if (distProb(gen) < mutationRate) {
            mutatedSeq[i] = allProcs[distProc(gen)].name;
        }
    }

    // Ajouter ces stratégies de mutation plus diversifiées
    // 10% chance de duplication d'une sous-séquence (générique, pas spécifique)
    if (distProb(gen) < 0.1 && mutatedSeq.size() > 5) {
        size_t startPos = gen() % (mutatedSeq.size() - 5);
        size_t length = 3 + gen() % 5; // Entre 3 et 7

        // Éviter de dépasser la taille
        if (startPos + length > mutatedSeq.size()) {
            length = mutatedSeq.size() - startPos;
        }

        vector<string> toRepeat(mutatedSeq.begin() + startPos, mutatedSeq.begin() + startPos + length);

        // Insérer la répétition à un endroit aléatoire
        size_t insertPos = gen() % (mutatedSeq.size() + 1);
        mutatedSeq.insert(mutatedSeq.begin() + insertPos, toRepeat.begin(), toRepeat.end());
    }

    // 10% chance d'inversion d'une sous-séquence
    if (distProb(gen) < 0.1 && mutatedSeq.size() > 5) {
        size_t startPos = gen() % (mutatedSeq.size() - 5);
        size_t length = 3 + gen() % 5; // Entre 3 et 7

        // Éviter de dépasser la taille
        if (startPos + length > mutatedSeq.size()) {
            length = mutatedSeq.size() - startPos;
        }

        std::reverse(mutatedSeq.begin() + startPos, mutatedSeq.begin() + startPos + length);
    }

    return Individual(mutatedSeq);
}

void ProcessManager::runGeneticAlgorithm() {
    // Ajouter au début de la méthode:
    auto startTime = std::chrono::high_resolution_clock::now();
    const int MAX_RUNTIME_SECONDS = 900; // 15 minutes maximum

    initializePopulation();
    evaluateFitness();

    // Track best fitness to potentially implement early stopping
    double bestFitnessSoFar = std::numeric_limits<double>::lowest();
    int generationsWithoutImprovement = 0;

    for (int generation = 0; generation < MAX_GENERATIONS; ++generation) {
        // Vérifier le temps écoulé au début de chaque génération
        auto currentTime = std::chrono::high_resolution_clock::now();
        auto elapsedSeconds = std::chrono::duration_cast<std::chrono::seconds>(currentTime - startTime).count();

        if (elapsedSeconds > MAX_RUNTIME_SECONDS) {
            std::cout << "Stopping after " << elapsedSeconds << " seconds (time limit reached)" << std::endl;
            break;
        }

        // Afficher progression toutes les 10 générations
        if (generation % 10 == 0) {
            std::cout << "Generation " << generation << ", best fitness: " << bestFitnessSoFar
                      << ", elapsed: " << elapsedSeconds << "s" << std::endl;
        }

        // Create the new population
        std::vector<Individual> newPopulation;

        // Elitism: Keep the best individuals
        std::vector<Individual> sortedPopulation = population_;
        std::sort(sortedPopulation.begin(), sortedPopulation.end(),
                  [](const Individual &a, const Individual &b) { return a.fitness > b.fitness; });

        size_t elitismCount = static_cast<size_t>(POPULATION_SIZE * ELITISM_RATE);
        for (size_t i = 0; i < elitismCount && i < sortedPopulation.size(); ++i) {
            newPopulation.push_back(sortedPopulation[i]);
        }

        // Tournament selection
        std::vector<Individual> parents = tournamentSelection(population_, POPULATION_SIZE - newPopulation.size());

        // Crossover and mutation
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dist(0.0, 1.0);

        while (newPopulation.size() < POPULATION_SIZE && parents.size() >= 2) {
            // Select two parents
            size_t idx1 = gen() % parents.size();
            size_t idx2;
            do {
                idx2 = gen() % parents.size();
            } while (idx2 == idx1 && parents.size() > 1);

            // Create children through crossover with probability CROSSOVER_RATE
            Individual child1 = parents[idx1];
            Individual child2 = parents[idx2];

            if (dist(gen) < CROSSOVER_RATE) {
                auto children = crossover(parents[idx1], parents[idx2]);
                child1 = children.first;
                child2 = children.second;
            }

            // Apply mutation
            child1 = mutate(child1, MUTATION_RATE);
            child2 = mutate(child2, MUTATION_RATE);

            // Add to new population
            if (newPopulation.size() < POPULATION_SIZE) {
                newPopulation.push_back(child1);
            }
            if (newPopulation.size() < POPULATION_SIZE) {
                newPopulation.push_back(child2);
            }
        }

        // Replace the old population
        population_ = newPopulation;

        // Evaluate fitness for new population
        evaluateFitness();

        // Si pas d'amélioration pendant longtemps, introduire des mutations destructrices
        // Pas de mutation forte avant au moins 50 générations
        if (generationsWithoutImprovement > 50) {
            // Appliquer une mutation plus agressive à 20% de la population (au lieu de 10%)
            int numToMutate = POPULATION_SIZE * 0.2;
            std::uniform_int_distribution<> distIndex(0, POPULATION_SIZE - 1);

            for (int i = 0; i < numToMutate; ++i) {
                int idx = distIndex(gen);
                // Doubler le taux de mutation normal pour ces individus
                population_[idx] = mutate(population_[idx], MUTATION_RATE * 2.0);
            }

            evaluateFitness(); // Ré-évaluer la population après les mutations
        }

        // Update best solution
        Individual currentBest = findBestIndividual(population_);
        if (currentBest.fitness > bestFitnessSoFar) {
            bestFitnessSoFar = currentBest.fitness;
            bestSolution_ = currentBest;
            generationsWithoutImprovement = 0;
        } else {
            generationsWithoutImprovement++;
        }

        // Dans la méthode runGeneticAlgorithm(), après avoir trouvé currentBest
        updateHallOfFame(currentBest);

        // Réinitialisation partielle si stagnation prolongée
        if (generationsWithoutImprovement % 30 == 0 && generationsWithoutImprovement > 0) {
            std::cout << "Applying stronger mutations instead of partial reset..." << std::endl;

            // Au lieu de réinitialiser, appliquer des mutations très fortes
            for (auto &individual : population_) {
                // 50% de chance d'appliquer une mutation forte
                if (dist(gen) < 0.5) {
                    individual = mutate(individual, 0.8); // Taux de mutation très élevé
                }
            }

            evaluateFitness();
        }

        // Early stopping if no improvement for 100 generations (au lieu de 30)
        if (generationsWithoutImprovement >= 100) {
            std::cout << "Stopping early after " << generationsWithoutImprovement << " generations without improvement"
                      << std::endl;
            break;
        }

        // Dans runGeneticAlgorithm(), à chaque 50 générations sans amélioration:

        if (generationsWithoutImprovement % 50 == 0 && generationsWithoutImprovement > 0) {
            std::cout << "Reintroducing diversity..." << std::endl;

            // Conserver les 20% meilleurs individus
            std::sort(population_.begin(), population_.end(),
                      [](const Individual &a, const Individual &b) { return a.fitness > b.fitness; });

            size_t keepCount = POPULATION_SIZE * 0.2;
            std::vector<Individual> elites(population_.begin(), population_.begin() + keepCount);

            // Réinitialiser avec de nouveaux individus
            initializePopulation();

            // Réintégrer l'élite
            for (size_t i = 0; i < elites.size(); ++i) {
                population_[i] = elites[i];
            }

            evaluateFitness();
        }
    }

    // Et à la fin de la méthode, avant generateOutput()
    // Utilisez le meilleur individu du hall of fame comme solution finale
    if (!hallOfFame_.empty()) {
        bestSolution_ = hallOfFame_[0];
    }

    generateOutput();
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

std::vector<Individual> ProcessManager::tournamentSelection(const std::vector<Individual> &population,
                                                            size_t numToSelect) {
    if (population.empty()) return {};

    std::vector<Individual> selected;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> dist(0, population.size() - 1);

    for (size_t i = 0; i < numToSelect; ++i) {
        Individual best = population[dist(gen)];

        // Select TOURNAMENT_SIZE random individuals and find the best
        for (int j = 1; j < TOURNAMENT_SIZE; ++j) {
            size_t idx = dist(gen);
            if (population[idx].fitness > best.fitness) {
                best = population[idx];
            }
        }

        selected.push_back(best);
    }

    return selected;
}

void ProcessManager::updateHallOfFame(const Individual &individual) {
    // Si le hall of fame n'est pas plein, ajoutez simplement l'individu
    if (hallOfFame_.size() < HALL_OF_FAME_SIZE) {
        hallOfFame_.push_back(individual);
    } else {
        // Trouver l'individu avec la fitness la plus basse dans le hall of fame
        auto minIt = std::min_element(hallOfFame_.begin(), hallOfFame_.end(),
                                      [](const Individual &a, const Individual &b) { return a.fitness < b.fitness; });

        // Si l'individu actuel est meilleur, remplacez le plus faible
        if (individual.fitness > minIt->fitness) {
            *minIt = individual;
        }
    }

    // Trier le hall of fame par fitness décroissante
    std::sort(hallOfFame_.begin(), hallOfFame_.end(),
              [](const Individual &a, const Individual &b) { return a.fitness > b.fitness; });
}
