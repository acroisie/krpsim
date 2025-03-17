#pragma once

#include "Config.hpp"
#include "Simulator.hpp"
#include <vector>
#include <string>
#include <functional>

class Optimizer {
public:
    struct Solution {
        std::vector<std::string> processSequence;
        double score;
    };

    Optimizer(const Config& config, int timeLimit);

    // Exécute l'algorithme d'optimisation et retourne la meilleure solution
    Solution optimize();

    // Afficher les résultats
    void printResults(const Solution& solution);

private:
    const Config& config_;
    int timeLimit_;
    Simulator simulator_;

    // Stratégies d'optimisation
    Solution greedyOptimize();
    Solution simulatedAnnealing(const Solution& initialSolution, int iterations);
    Solution hillClimbing(const Solution& initialSolution, int iterations);
    Solution randomSearch(int numTrials);

    // Opérations sur les solutions
    Solution randomSolution();
    Solution neighborSolution(const Solution& solution);
    
    // Fonctions utilitaires
    double evaluate(const Solution& solution);
    static bool compareSolutions(const Solution& a, const Solution& b);
    
    // Analyse et stratégies
    std::vector<std::string> analyzeProcessDependencies();
    std::vector<std::string> prioritizeByResourceEfficiency();
};