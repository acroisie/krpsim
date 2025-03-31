#include "ProcessManager.hpp"
#include <iostream>
#include <limits>

using namespace std;

ProcessManager::ProcessManager(const Config &config, int delayLimit)
    : config_(config), simulator_(config, delayLimit), geneticAlgorithm_(config, simulator_) {
    // Initialiser les stocks actuels
    for (const auto &stock : config_.getStocks()) {
        currentStocks_[stock.name] = stock.quantity;
    }
}

void ProcessManager::run() {
    cout << "Nice file! " << config_.getProcesses().size() << " processes, " << config_.getStocks().size()
         << " stocks, " << config_.getOptimizeGoal().size() << " to optimize" << endl;

    cout << "Evaluating .................. " << flush;

    // Exécuter l'algorithme génétique
    Individual bestSolution = geneticAlgorithm_.runGeneticAlgorithm();

    cout << "done." << endl;

    // Générer la sortie
    generateOutput(bestSolution);
}

void ProcessManager::generateOutput(const Individual &bestSolution) {
    // Reset state for final simulation
    currentStocks_.clear();
    for (const auto &stock : config_.getStocks()) {
        currentStocks_[stock.name] = stock.quantity;
    }
    executionLogs_.clear();

    // Calculer le fitness pour obtenir les logs d'exécution
    simulator_.calculateFitness(bestSolution.processSequence, executionLogs_);

    // Afficher les résultats
    cout << "Main walk" << endl;
    for (auto &log : executionLogs_) {
        cout << log.first << ":" << log.second << endl;
    }

    // Simuler la séquence pour obtenir les stocks finaux
    map<string, int> finalStocks;
    int executedCount;
    simulator_.simulateSequence(bestSolution.processSequence, finalStocks, executedCount);

    // Déterminer le cycle final (temps max)
    int finalCycle = 0;
    if (!executionLogs_.empty()) {
        // Trouver le dernier processus exécuté
        int lastStartCycle = executionLogs_.back().first;
        string lastProcessName = executionLogs_.back().second;

        // Trouver le processus et sa durée
        for (const auto &proc : config_.getProcesses()) {
            if (proc.name == lastProcessName) {
                finalCycle = lastStartCycle + proc.nbCycle;
                break;
            }
        }
    }

    cout << "no more process doable at time " << finalCycle << endl;

    cout << "Stock :" << endl;
    for (auto &stock : finalStocks) {
        cout << stock.first << " => " << stock.second << endl;
    }
}