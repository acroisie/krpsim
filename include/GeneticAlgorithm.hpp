#pragma once
#include "Config.hpp"
#include "Individual.hpp"
#include "Simulator.hpp"
#include <vector>
#include <string>

class GeneticAlgorithm {
public:
    GeneticAlgorithm(const Config& config, Simulator& simulator, int populationSize = 50);
    
    // Initialise la population avec des individus valides
    void initializePopulation();
    
    // Évalue le fitness de tous les individus
    void evaluateFitness();
    
    // Sélectionne des individus pour la reproduction
    std::vector<Individual> selectParents(const std::vector<Individual>& population);
    
    // Effectue un croisement entre deux parents
    std::pair<Individual, Individual> crossover(const Individual& parent1, const Individual& parent2);
    
    // Effectue une mutation sur un individu
    Individual mutate(const Individual& individual, double mutationRate);
    
    // Trouve le meilleur individu dans une population
    Individual findBestIndividual(const std::vector<Individual>& population);
    
    // Exécute l'algorithme génétique
    Individual runGeneticAlgorithm(int numGenerations = 50);
    
    // Accesseur pour les noms de tous les processus
    const std::vector<std::string>& getAllProcessNames() const;

private:
    const Config& config_;
    Simulator& simulator_;
    std::vector<Individual> population_;
    int populationSize_;
    std::vector<std::string> allProcessNames_;
};