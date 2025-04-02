#include "ProcessManager.hpp"
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <limits> // Pour std::numeric_limits
#include <queue>  // Pour std::priority_queue dans la re-simulation
#include <algorithm> // Pour std::min

// Constructeur
ProcessManager::ProcessManager(const Config& config, int delayLimit)
    : config_(config),
      simulator_(config, delayLimit), // Initialise le simulateur
      // Initialise l'AG avec des valeurs par défaut (peuvent être ajustées)
      geneticAlgorithm_(config, simulator_, 100, 0.05, 0.7, 2, 50, 150),
      delayLimit_(delayLimit)
{}

// Lance l'exécution
void ProcessManager::run() {
    // Affichage initial
    std::cout << "Nice file! "
              << config_.getProcesses().size() << " processes, "
              << config_.getStocks().size() << " initial stocks, "
              << config_.getOptimizeGoal().size() << " optimization goal(s)" << std::endl;

    std::cout << "Evaluating using Genetic Algorithm..." << std::endl;

    // --- Exécution de l'Algorithme Génétique ---
    // Définir le nombre de générations (peut être un paramètre)
    int numberOfGenerations = 100; // Exemple, à ajuster
    Individual bestSolution = geneticAlgorithm_.runEvolution(numberOfGenerations);

    std::cout << "Optimization complete." << std::endl;

    // --- Génération de la Sortie ---
    generateOutput(bestSolution);
}

// Génère la sortie finale
void ProcessManager::generateOutput(const Individual& bestSolution) {
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Best solution found:" << std::endl;

    // Ré-exécuter la simulation une dernière fois avec la meilleure séquence
    // pour obtenir les logs d'exécution finaux.
    std::vector<std::pair<int, std::string>> finalExecutionLogs;

    // Recalculer le fitness pour obtenir les logs propres de la meilleure solution.
    double finalFitness = simulator_.calculateFitnessAndLogs(
        bestSolution.processExecutionAttemptSequence,
        finalExecutionLogs);

     std::cout << "Final Fitness Score: " << finalFitness << std::endl;
     std::cout << "Execution Log (Format: <cycle>:<process_name>):" << std::endl;

    // Afficher les logs d'exécution et calculer le temps final
    int lastCycle = 0;
    if (finalExecutionLogs.empty()) {
        std::cout << "(No processes executed)" << std::endl;
    } else {
        for (const auto& log : finalExecutionLogs) {
            std::cout << log.first << ":" << log.second << std::endl;
             const Process* currentProc = nullptr;
              for (const auto& p : config_.getProcesses()) {
                  if (p.name == log.second) {
                      currentProc = &p;
                      break;
                  }
              }
             if (currentProc) {
                 int endTime = log.first + currentProc->nbCycle;
                 if (endTime > lastCycle) {
                     lastCycle = endTime;
                 }
             }
        }
    }
     lastCycle = std::min(lastCycle, delayLimit_);


    std::cout << "----------------------------------------" << std::endl;
    // Déterminer le message final basé sur lastCycle
    if (finalExecutionLogs.empty()) {
         std::cout << "No process could be executed within the time limit (" << delayLimit_ << ")." << std::endl;
    } else if (lastCycle < delayLimit_) {
         std::cout << "Simulation ended at cycle " << lastCycle << "." << std::endl;
    } else {
         std::cout << "Simulation reached time limit at cycle " << delayLimit_ << "." << std::endl;
    }


    // --- Re-simulation pour obtenir les stocks finaux (logique simplifiée) ---
    // *** TEMPORAIREMENT COMMENTÉ POUR ISOLER L'ERREUR ASan ***
    /*
    std::cout << "Calculating final stocks (simplified re-simulation)..." << std::endl;
    std::map<std::string, int> stocksAtEnd;
    for (const auto& stock : config_.getStocks()) {
        stocksAtEnd[stock.name] = stock.quantity;
    }
    int simTime = 0;
    std::priority_queue<Simulator::RunningProcessInfo, std::vector<Simulator::RunningProcessInfo>, std::greater<Simulator::RunningProcessInfo>> runningQueue;

     std::vector<std::pair<int, const Process*>> loggedProcessPtrs;
     loggedProcessPtrs.reserve(finalExecutionLogs.size());
     for(const auto& log : finalExecutionLogs) {
         const Process* p_ptr = nullptr;
         for(const auto& p : config_.getProcesses()) {
             if(p.name == log.second) {
                 p_ptr = &p;
                 break;
             }
         }
         if (p_ptr) {
            loggedProcessPtrs.push_back({log.first, p_ptr});
         } else {
             std::cerr << "Warning: Logged process '" << log.second << "' not found in config during final stock calculation." << std::endl;
         }
     }

     size_t logIdx = 0;
     while(simTime < delayLimit_) {
         while(!runningQueue.empty() && runningQueue.top().completionTime <= simTime) {
             auto finished = runningQueue.top();
             runningQueue.pop();
             if(finished.processDetails) {
                 for(const auto& out : finished.processDetails->outputs) {
                     stocksAtEnd[out.first] += out.second;
                 }
             }
         }

         while(logIdx < loggedProcessPtrs.size() && loggedProcessPtrs[logIdx].first == simTime) {
              const Process* procToStart = loggedProcessPtrs[logIdx].second;
              for(const auto& in : procToStart->inputs) {
                   if (stocksAtEnd.count(in.first)) {
                       stocksAtEnd[in.first] -= in.second;
                   } else {
                        stocksAtEnd[in.first] = -in.second;
                        std::cerr << "Warning: Consuming resource '" << in.first
                                  << "' which was not initially present during final stock calculation." << std::endl;
                   }
              }
              runningQueue.push({procToStart, simTime + procToStart->nbCycle});
              logIdx++;
         }

         int nextEventTime = delayLimit_;
         if (!runningQueue.empty()) {
             nextEventTime = std::min(nextEventTime, runningQueue.top().completionTime);
         }
         if (logIdx < loggedProcessPtrs.size()) {
             nextEventTime = std::min(nextEventTime, loggedProcessPtrs[logIdx].first);
         }

         if (nextEventTime <= simTime || nextEventTime >= delayLimit_) {
             break;
         }
         simTime = nextEventTime;
     }


    std::cout << "Final Stocks:" << std::endl;
    std::map<std::string, bool> stockPrinted;
     for (const auto& stockPair : stocksAtEnd) {
        std::cout << "  " << stockPair.first << ": " << stockPair.second << std::endl;
        stockPrinted[stockPair.first] = true;
    }
     for (const auto& proc : config_.getProcesses()) {
         for (const auto& input : proc.inputs) {
             if (stockPrinted.find(input.first) == stockPrinted.end()) {
                 std::cout << "  " << input.first << ": 0" << std::endl;
                 stockPrinted[input.first] = true;
             }
         }
         for (const auto& output : proc.outputs) {
              if (stockPrinted.find(output.first) == stockPrinted.end()) {
                 std::cout << "  " << output.first << ": 0" << std::endl;
                 stockPrinted[output.first] = true;
             }
         }
     }
      for (const auto& initialStock : config_.getStocks()) {
           if (stockPrinted.find(initialStock.name) == stockPrinted.end()) {
                std::cout << "  " << initialStock.name << ": " << initialStock.quantity << " (initial, unused)" << std::endl;
                stockPrinted[initialStock.name] = true;
           }
      }
    */
     // Message temporaire pendant que la re-simulation est commentée
     std::cout << "Final Stocks: (Calculation temporarily disabled for debugging)" << std::endl;


    std::cout << "----------------------------------------" << std::endl;
}
