#include "ProcessManager.hpp"
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <limits>
#include <algorithm> // Pour std::max utilisé dans le calcul du temps final affiché
#include <set>       // Pour afficher tous les stocks pertinents

// Constructeur (inchangé)
ProcessManager::ProcessManager(const Config& config, int delayLimit)
    : config_(config),
      simulator_(config, delayLimit),
      geneticAlgorithm_(config, simulator_, 100, 0.05, 0.7, 2, 50, 150),
      delayLimit_(delayLimit)
{}

// Lance l'exécution (inchangé)
void ProcessManager::run() {
    std::cout << "Nice file! "
              << config_.getProcesses().size() << " processes, "
              << config_.getStocks().size() << " initial stocks, "
              << config_.getOptimizeGoal().size() << " optimization goal(s)" << std::endl;
    std::cout << "Evaluating using Genetic Algorithm..." << std::endl;
    int numberOfGenerations = 100;
    Individual bestSolution = geneticAlgorithm_.runEvolution(numberOfGenerations);
    std::cout << "Optimization complete." << std::endl;
    generateOutput(bestSolution);
}

// Génère la sortie finale en utilisant SimulationResult
void ProcessManager::generateOutput(const Individual& bestSolution) {
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Best solution found:" << std::endl;

    // *** Exécuter la simulation UNE SEULE FOIS pour obtenir tous les résultats ***
    Simulator::SimulationResult finalResult = simulator_.calculateFitnessAndLogs(
        bestSolution.processExecutionAttemptSequence);

    // Utiliser les résultats retournés
    double finalFitness = finalResult.fitness;
    const auto& finalExecutionLogs = finalResult.logs;
    const auto& stocksAtEnd = finalResult.finalStocks; // Stocks finaux corrects
    int finalCycle = finalResult.finalCycle;           // Temps final correct
    bool reachedTimeLimit = finalResult.reachedTimeLimit; // Indicateur de limite atteinte

     std::cout << "Final Fitness Score: " << finalFitness << std::endl;
     std::cout << "Execution Log (Format: <cycle>:<process_name>):" << std::endl;

    // Afficher les logs d'exécution (inchangé)
    if (finalExecutionLogs.empty()) {
        std::cout << "(No processes executed)" << std::endl;
    } else {
        for (const auto& log : finalExecutionLogs) {
            std::cout << log.first << ":" << log.second << std::endl;
        }
    }

    std::cout << "----------------------------------------" << std::endl;

    // Afficher le message de fin basé sur le résultat de la simulation
    // Utiliser le format de l'exemple KrpSim si possible
    if (finalExecutionLogs.empty()) {
         std::cout << "No process could be executed within the time limit (" << delayLimit_ << ")." << std::endl;
    } else if (!reachedTimeLimit && finalCycle < delayLimit_) {
         // Si on n'a pas atteint la limite et que le dernier événement est avant,
         // on peut supposer qu'aucun autre processus n'était possible.
         std::cout << "no more process doable at time " << finalCycle << std::endl;
    }
     else {
         // Si on a atteint la limite ou si le dernier événement se termine pile à la limite
          std::cout << "Simulation reached time limit at cycle " << delayLimit_ << "." << std::endl;
          // Ou on pourrait aussi utiliser le message "no more process doable..." s'il n'y a plus rien dans la queue à la fin?
          // Gardons simple pour l'instant.
    }


    // --- Affichage des Stocks Finaux (Maintenant Correct) ---
    std::cout << "Stock:" << std::endl; // Format requis par le sujet

    // Pour afficher tous les stocks pertinents (initiaux + produits), même à zéro.
    std::set<std::string> relevantStocks;
    // Ajouter les stocks initiaux
    for(const auto& stock : config_.getStocks()) {
        relevantStocks.insert(stock.name);
    }
    // Ajouter les stocks mentionnés dans les processus (inputs/outputs)
    for(const auto& proc : config_.getProcesses()) {
        for(const auto& input : proc.inputs) relevantStocks.insert(input.first);
        for(const auto& output : proc.outputs) relevantStocks.insert(output.first);
    }

    // Afficher chaque stock pertinent avec sa valeur finale
    for(const std::string& stockName : relevantStocks) {
        int finalQuantity = 0;
        auto it = stocksAtEnd.find(stockName);
        if (it != stocksAtEnd.end()) {
            finalQuantity = it->second;
        }
        // Afficher au format "nom => quantité"
        std::cout << "  " << stockName << " => " << finalQuantity << std::endl;
    }

    std::cout << "----------------------------------------" << std::endl;
}
