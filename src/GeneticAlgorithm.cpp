#include "GeneticAlgorithm.hpp"
#include <algorithm>
#include <iostream>
#include <limits>
#include <random>
#include <cmath>
#include <unordered_set>
#include <queue>

using namespace std;

const size_t MAX_SEQUENCE_LENGTH = 500;


GeneticAlgorithm::GeneticAlgorithm(const Config& config, Simulator& simulator, int populationSize)
    : config_(config), simulator_(simulator), populationSize_(populationSize) {
    
    // Collect all process names
    for (const auto& proc : config_.getProcesses()) {
        allProcessNames_.push_back(proc.name);
    }
}

const vector<string>& GeneticAlgorithm::getAllProcessNames() const {
    return allProcessNames_;
}

void GeneticAlgorithm::initializePopulation() {
    population_.resize(populationSize_);
    
    random_device rd;
    mt19937 gen(rd());
    
    // Créer différents types d'individus pour augmenter la diversité
    for (int i = 0; i < populationSize_; ++i) {
        vector<string> sequence;
        
        // Varier la longueur des séquences pour explorer différentes options
        int seqLen;
        if (i % 3 == 0) {
            // Un tiers des individus ont des séquences courtes
            seqLen = 30 + (gen() % 20); // 30-50 instructions
        } else if (i % 3 == 1) {
            // Un tiers ont des séquences moyennes
            seqLen = 50 + (gen() % 30); // 50-80 instructions
        } else {
            // Un tiers ont des séquences longues
            seqLen = 80 + (gen() % 40); // 80-120 instructions
        }
        
        sequence.reserve(seqLen);
        
        // Générer une séquence aléatoire
        for (int j = 0; j < seqLen; ++j) {
            uniform_int_distribution<> procDist(0, static_cast<int>(allProcessNames_.size()) - 1);
            sequence.push_back(allProcessNames_[procDist(gen)]);
        }
        
        population_[i] = Individual(sequence);
    }
}

void GeneticAlgorithm::evaluateFitness() {
    for (auto &indiv : population_) {
        vector<pair<int, string>> logs;
        indiv.fitness = simulator_.calculateFitness(indiv.processSequence, logs);
    }
}

vector<Individual> GeneticAlgorithm::selectParents(const vector<Individual> &population) {
    vector<Individual> parents;
    if (population.empty()) {
        return parents;
    }

    double totalFitness = 0.0;
    double minFitness = numeric_limits<double>::max();

    for (const auto &indiv : population) {
        if (indiv.fitness < minFitness) {
            minFitness = indiv.fitness;
        }
    }

    // On shift si jamais toutes les fitness sont négatives
    vector<double> adjustedFitness;
    adjustedFitness.reserve(population.size());
    
    // Shift pour avoir des valeurs positives
    double shift = std::abs(minFitness) + 1.0;
    totalFitness = 0.0;
    for (const auto &ind : population) {
        double adj = ind.fitness + shift;
        adjustedFitness.push_back(adj);
        totalFitness += adj;
    }

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dist(0.0, totalFitness);

    int nbParents = min(populationSize_, static_cast<int>(population.size()));
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
        // Fallback si rien n'est sélectionné
        if (static_cast<int>(parents.size()) <= i && !population.empty()) {
            parents.push_back(population[0]);
        }
    }
    return parents;
}

pair<Individual, Individual> GeneticAlgorithm::crossover(const Individual &parent1, const Individual &parent2) {
    if (parent1.processSequence.empty() || parent2.processSequence.empty()) {
        return make_pair(parent1, parent2);
    }
    
    random_device rd;
    mt19937 gen(rd());
    
    size_t minSize = min(parent1.processSequence.size(), parent2.processSequence.size());
    if (minSize < 2) {
        return make_pair(parent1, parent2);
    }

    // Utiliser un double point de croisement pour plus de variété génétique
    uniform_int_distribution<> pointDist1(1, static_cast<int>(minSize) / 3);
    uniform_int_distribution<> pointDist2(static_cast<int>(minSize) / 3 + 1, static_cast<int>(minSize) - 1);
    
    size_t crossoverPoint1 = pointDist1(gen);
    size_t crossoverPoint2 = pointDist2(gen);
    
    vector<string> child1, child2;
    
    // Premier enfant: début parent1, milieu parent2, fin parent1
    child1.insert(child1.end(), 
                 parent1.processSequence.begin(), 
                 parent1.processSequence.begin() + crossoverPoint1);
    
    if (crossoverPoint1 < parent2.processSequence.size() && crossoverPoint2 < parent2.processSequence.size()) {
        child1.insert(child1.end(), 
                    parent2.processSequence.begin() + crossoverPoint1, 
                    parent2.processSequence.begin() + crossoverPoint2);
    }
    
    if (crossoverPoint2 < parent1.processSequence.size()) {
        child1.insert(child1.end(), 
                    parent1.processSequence.begin() + crossoverPoint2, 
                    parent1.processSequence.end());
    }
    
    // Deuxième enfant: début parent2, milieu parent1, fin parent2
    child2.insert(child2.end(), 
                 parent2.processSequence.begin(), 
                 parent2.processSequence.begin() + crossoverPoint1);
    
    if (crossoverPoint1 < parent1.processSequence.size() && crossoverPoint2 < parent1.processSequence.size()) {
        child2.insert(child2.end(), 
                    parent1.processSequence.begin() + crossoverPoint1, 
                    parent1.processSequence.begin() + crossoverPoint2);
    }
    
    if (crossoverPoint2 < parent2.processSequence.size()) {
        child2.insert(child2.end(), 
                    parent2.processSequence.begin() + crossoverPoint2, 
                    parent2.processSequence.end());
    }

    return make_pair(Individual(child1), Individual(child2));
}

Individual GeneticAlgorithm::mutate(const Individual &individual, double mutationRate) {
    auto allProcs = config_.getProcesses();
    if (allProcs.empty() || individual.processSequence.empty()) {
        return individual;
    }

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> distProb(0.0, 1.0);
    
    vector<string> mutatedSeq = individual.processSequence;
    
    // Mutation ponctuelle - remplacer des processus aléatoires
    for (size_t i = 0; i < mutatedSeq.size(); ++i) {
        if (distProb(gen) < mutationRate) {
            uniform_int_distribution<> procDist(0, static_cast<int>(allProcessNames_.size()) - 1);
            mutatedSeq[i] = allProcessNames_[procDist(gen)];
        }
    }
    
    // Mutation structurelle - ajouter ou supprimer des segments
    if (distProb(gen) < 0.3) {
        // 30% de chance d'avoir une mutation structurelle
        
        if (distProb(gen) < 0.5 && mutatedSeq.size() > 10) {
            // Supprimer un segment
            uniform_int_distribution<> startDist(0, static_cast<int>(mutatedSeq.size()) - 10);
            uniform_int_distribution<> lengthDist(1, 10);
            
            int start = startDist(gen);
            int length = min(lengthDist(gen), static_cast<int>(mutatedSeq.size()) - start);
            
            mutatedSeq.erase(mutatedSeq.begin() + start, mutatedSeq.begin() + start + length);
        } else {
            // Ajouter un segment
            uniform_int_distribution<> posDist(0, static_cast<int>(mutatedSeq.size()));
            uniform_int_distribution<> lengthDist(1, 10);
            
            int pos = posDist(gen);
            int length = lengthDist(gen);
            
            for (int i = 0; i < length; ++i) {
                uniform_int_distribution<> procDist(0, static_cast<int>(allProcessNames_.size()) - 1);
                mutatedSeq.insert(mutatedSeq.begin() + pos, allProcessNames_[procDist(gen)]);
            }
        }
    }
    
    return Individual(mutatedSeq);
}

Individual GeneticAlgorithm::findBestIndividual(const vector<Individual> &population) {
    if (population.empty()) {
        return Individual();
    }
    
    Individual best = population.front();
    for (const auto &ind : population) {
        if (ind.fitness > best.fitness) {
            best = ind;
        }
    }
    return best;
}

Individual GeneticAlgorithm::runGeneticAlgorithm(int numGenerations) {
    cout << "Starting genetic algorithm optimization..." << endl;
    initializePopulation();
    
    evaluateFitness();
    
    // Track the best solution found
    Individual bestSolution;
    if (!population_.empty()) {
        bestSolution = findBestIndividual(population_);
    }

    const double MUTATION_RATE = 0.15;
    
    // Run for specified number of generations
    for (int generation = 0; generation < numGenerations; ++generation) {
        vector<Individual> parents = selectParents(population_);
        vector<Individual> newPopulation;

        // Elitism - keep the best solution and some more top solutions
        if (!population_.empty()) {
            // Trier la population par fitness
            vector<Individual> sortedPopulation = population_;
            sort(sortedPopulation.begin(), sortedPopulation.end(), 
                 [](const Individual& a, const Individual& b) { 
                     return a.fitness > b.fitness; 
                 });
                 
            // Garder les 3 meilleurs individus (élitisme)
            int eliteCount = min(3, static_cast<int>(sortedPopulation.size()));
            for (int i = 0; i < eliteCount; i++) {
                newPopulation.push_back(sortedPopulation[i]);
                
                // Mettre à jour la meilleure solution
                if (i == 0 && sortedPopulation[0].fitness > bestSolution.fitness) {
                    bestSolution = sortedPopulation[0];
                }
            }
        }

        // Create new population through crossover and mutation
        while (newPopulation.size() < population_.size() && parents.size() >= 2) {
            // Select two parents randomly
            random_device rd;
            mt19937 gen(rd());
            uniform_int_distribution<> dist(0, static_cast<int>(parents.size()) - 1);
            
            size_t idx1 = dist(gen);
            size_t idx2 = dist(gen);
            while (idx2 == idx1 && parents.size() > 1) {
                idx2 = dist(gen);
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

        // If we don't have enough individuals, add random ones
        while (newPopulation.size() < population_.size()) {
            random_device rd;
            mt19937 gen(rd());
            
            vector<string> randomSeq;
            uniform_int_distribution<> lenDist(30, 120);
            int seqLen = lenDist(gen);
            
            for (int j = 0; j < seqLen; j++) {
                uniform_int_distribution<> procDist(0, static_cast<int>(allProcessNames_.size()) - 1);
                randomSeq.push_back(allProcessNames_[procDist(gen)]);
            }
            
            newPopulation.push_back(Individual(randomSeq));
        }

        // Replace population
        if (!newPopulation.empty()) {
            population_ = newPopulation;
        }
        
        // Evaluate fitness of new population
        evaluateFitness();
        
        // Track the best solution in this generation
        auto currentBest = findBestIndividual(population_);
        
        // Update best solution if improvement found
        if (currentBest.fitness > bestSolution.fitness) {
            bestSolution = currentBest;
            cout << "Generation " << generation + 1 << ": New best fitness = " << bestSolution.fitness << endl;
        } else {
            cout << "Generation " << generation + 1 << ": No improvement (best = " << bestSolution.fitness << ")" << endl;
        }
    }

    cout << "Genetic algorithm completed." << endl;
    return bestSolution;
}