#include "ProcessManager.hpp"
#include <iostream>

ProcessManager::ProcessManager(const Config &config, int delayLimit)
    : config_(config), simulator_(config), optimizer_(config, delayLimit) {
}

void ProcessManager::runGeneticAlgorithm() {
    std::cout << "Starting optimization..." << std::endl;
    
    // Exécuter l'optimisation
    Optimizer::Solution bestSolution = optimizer_.optimize();
    
    // Afficher les résultats
    optimizer_.printResults(bestSolution);
}