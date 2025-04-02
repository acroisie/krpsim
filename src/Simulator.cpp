#include "Simulator.hpp"
#include <iostream>
#include <limits>
#include <queue>
#include <map>
#include <vector>
#include <string>
#include <set>
#include <algorithm> // Pour std::min, std::max
#include <cmath>     // Pour std::isfinite
#include <stdexcept> // Pour std::runtime_error

// Constructeur (inchangé)
Simulator::Simulator(const Config& config, int timeLimit)
    : systemConfig_(config), simulationTimeLimit_(timeLimit)
{
    processLookupMap_.clear();
    for (const auto& proc : systemConfig_.getProcesses()) {
        if (processLookupMap_.count(proc.name)) {
             throw std::runtime_error("Simulator Error: Duplicate process name found in config: " + proc.name);
        }
        processLookupMap_[proc.name] = &proc;
    }
     if (systemConfig_.getProcesses().empty()) {
          std::cerr << "Warning: Simulator initialized with no processes defined in the configuration." << std::endl;
     }
     // Debug: Vérifier si la map est peuplée
     // std::cout << "Debug: processLookupMap_ contains " << processLookupMap_.size() << " entries." << std::endl;
}

// Helper : Trouve un processus par son nom (inchangé)
const Process* Simulator::findProcessByName(const std::string& name) const {
    auto it = processLookupMap_.find(name);
    if (it != processLookupMap_.end()) { return it->second; }
    return nullptr;
}

// Helper : Vérifie si un processus peut démarrer (inchangé)
bool Simulator::canProcessStart(const Process* process, const std::map<std::string, int>& currentStocks) const {
    if (!process) { return false; }
    for (const auto& input : process->inputs) {
        auto it = currentStocks.find(input.first);
        if (it == currentStocks.end() || it->second < input.second) { return false; }
    }
    return true;
}

// Fonction CLÉ : Calcule le fitness et génère les logs d'exécution
Simulator::SimulationResult Simulator::calculateFitnessAndLogs(
    const std::vector<std::string>& sequenceToAttempt) {

    SimulationResult result;
    result.finalStocks.clear();
    for (const auto& stock : systemConfig_.getStocks()) {
        result.finalStocks[stock.name] = stock.quantity;
    }
    result.logs.clear();
    int currentTime = 0;
    int executedProcessCount = 0;
    size_t sequenceIndex = 0;

    std::priority_queue<RunningProcessInfo, std::vector<RunningProcessInfo>, std::greater<RunningProcessInfo>> runningProcessesQueue;
    bool opportunisticMode = false;

    while (currentTime <= simulationTimeLimit_) {

        // --- A. Compléter les processus ---
        bool resourcesUpdated = false;
        while (!runningProcessesQueue.empty() && runningProcessesQueue.top().completionTime <= currentTime) {
            RunningProcessInfo finishedProcess = runningProcessesQueue.top();
            runningProcessesQueue.pop();
            if (finishedProcess.completionTime <= simulationTimeLimit_ && finishedProcess.processDetails) {
                 for (const auto& output : finishedProcess.processDetails->outputs) {
                    result.finalStocks[output.first] += output.second;
                    resourcesUpdated = true;
                 }
            }
        }
        if (currentTime >= simulationTimeLimit_) {
             result.reachedTimeLimit = true;
             break;
        }

        // --- B. Essayer de lancer de nouveaux processus ---
        bool processStartedThisCycle = false;
        if (currentTime < simulationTimeLimit_) {
            bool tryStartingMore = true;
            while (tryStartingMore) {
                tryStartingMore = false;
                const Process* processToTry = nullptr;
                std::string processNameToLog;
                bool attemptedFromSequence = false;

                // Choisir quel processus essayer:
                if (!opportunisticMode && sequenceIndex < sequenceToAttempt.size()) {
                    processNameToLog = sequenceToAttempt[sequenceIndex];
                    processToTry = findProcessByName(processNameToLog);
                    attemptedFromSequence = true;
                    // NE PAS AVANCER sequenceIndex ici
                } else {
                    opportunisticMode = true; // Activer ou rester en mode opportuniste
                    const Process* bestOpportunisticChoice = nullptr;
                    for(const auto& pair : processLookupMap_) {
                        const Process* potentialProcess = pair.second;
                        if (canProcessStart(potentialProcess, result.finalStocks)) {
                            bestOpportunisticChoice = potentialProcess;
                            processNameToLog = potentialProcess->name;
                            break; // Stratégie simple: prendre le premier possible
                        }
                    }
                    processToTry = bestOpportunisticChoice;
                    if (!processToTry) break; // Rien à faire en opportuniste pour ce tour
                }

                // Si on a trouvé un processus à essayer
                if (processToTry) {
                    if (canProcessStart(processToTry, result.finalStocks) && (currentTime + processToTry->nbCycle <= simulationTimeLimit_)) {
                        // Lancer le processus
                        for (const auto& input : processToTry->inputs) {
                            result.finalStocks[input.first] -= input.second;
                        }
                        runningProcessesQueue.push({processToTry, currentTime + processToTry->nbCycle});
                        result.logs.push_back({currentTime, processNameToLog});
                        executedProcessCount++;
                        processStartedThisCycle = true;
                        tryStartingMore = true; // Essayer d'en lancer un autre

                        // *** CORRECTION : Avancer sequenceIndex SEULEMENT si succès ET depuis séquence ***
                        if (attemptedFromSequence) {
                            sequenceIndex++;
                        }
                    } else {
                        // Échec du lancement (ressources ou temps)
                        // *** CORRECTION : Ne PAS avancer sequenceIndex en cas d'échec ***
                        // if (attemptedFromSequence) {
                        //     sequenceIndex++; // <- Suppression de cette ligne
                        // }
                        // Si on a échoué à lancer depuis la séquence, on essaiera à nouveau au prochain cycle
                        // ou le mode opportuniste prendra le relais si la séquence finit.
                        // Si on est en mode opportuniste et qu'on échoue (ne devrait pas arriver), on sort.
                        if(!attemptedFromSequence) tryStartingMore = false;
                    }
                } else {
                    // processToTry est nullptr (nom invalide dans séquence ou fin de séquence)
                    if (attemptedFromSequence) {
                        // Nom invalide dans la séquence, il faut passer au suivant
                         sequenceIndex++; // <<< On avance ici pour ignorer l'étape invalide
                    }
                    // Vérifier si on passe en mode opportuniste
                    if (!opportunisticMode && sequenceIndex >= sequenceToAttempt.size()) {
                        opportunisticMode = true;
                        tryStartingMore = true; // Relancer la boucle interne pour tester l'opportuniste
                    }
                }
            } // Fin de la boucle interne while(tryStartingMore)
        } // Fin if (currentTime < simulationTimeLimit_)


        // --- C. Avancer le Temps ---
        if (processStartedThisCycle || resourcesUpdated) {
             if (currentTime >= simulationTimeLimit_) { result.reachedTimeLimit = true; break; }
             continue;
        } else {
            if (runningProcessesQueue.empty()) {
                break;
            } else {
                int nextCompletionTime = runningProcessesQueue.top().completionTime;
                if (nextCompletionTime > simulationTimeLimit_) {
                     currentTime = simulationTimeLimit_;
                     result.reachedTimeLimit = true;
                } else if (nextCompletionTime <= currentTime) {
                     // Sécurité anti-blocage
                     currentTime++;
                } else {
                    currentTime = nextCompletionTime;
                }
            }
        }
         if (currentTime > simulationTimeLimit_) {
              currentTime = simulationTimeLimit_;
              result.reachedTimeLimit = true;
              break;
         }
    } // Fin while (currentTime <= simulationTimeLimit_)

    // --- D. Dernière passe de complétion ---
     while (!runningProcessesQueue.empty() && runningProcessesQueue.top().completionTime <= simulationTimeLimit_) {
         RunningProcessInfo finishedProcess = runningProcessesQueue.top();
         runningProcessesQueue.pop();
         if (finishedProcess.processDetails) {
              for (const auto& output : finishedProcess.processDetails->outputs) {
                 result.finalStocks[output.first] += output.second;
              }
         }
     }

    // --- 3. Calcul Final du Fitness ---
    // (Le reste de la fonction est inchangé par rapport à v7)
    result.fitness = 0.0; // Renommé pour clarté
    const auto& optimizationGoals = systemConfig_.getOptimizeGoal();

    if (executedProcessCount == 0) {
        result.fitness = -1.0e18;
        result.finalCycle = currentTime;
        return result;
    }

    bool timeIsGoal = false;
    if (!result.logs.empty()) {
        int max_end_time = 0;
        for(const auto& log : result.logs) {
            const Process* proc = findProcessByName(log.second);
            if (proc) {
                max_end_time = std::max(max_end_time, log.first + proc->nbCycle);
            }
        }
        result.finalCycle = std::min(max_end_time, simulationTimeLimit_);
    } else {
         result.finalCycle = 0;
    }

    for (const std::string& goal : optimizationGoals) {
        if (goal == "time") {
            timeIsGoal = true;
            result.fitness += 10000.0 / (1.0 + result.finalCycle);
        } else {
            double stockValue = 0;
            if (result.finalStocks.count(goal)) {
                stockValue = static_cast<double>(result.finalStocks.at(goal));
            }
            double weight = 1.0;
            if (goal == "armoire") { weight = 100.0; }
            if (goal == "euro") { weight = 1.0; }
            if (goal == "boite") { weight = 5.0; }
            if (goal == "tarte_pomme" || goal == "tarte_citron" || goal == "flan") { weight = 0.1; }
            result.fitness += stockValue * weight;
        }
    }

    result.fitness += executedProcessCount * (simulationTimeLimit_ > 0 ? (50.0 / simulationTimeLimit_) : 0.01);

    if (!timeIsGoal && currentTime < simulationTimeLimit_ && runningProcessesQueue.empty() && opportunisticMode) {
         result.fitness *= 0.5;
    }

     if (!std::isfinite(result.fitness)) {
          std::cerr << "Warning: Calculated fitness is not finite (" << result.fitness << "). Resetting to large negative." << std::endl;
          result.fitness = -1.0e17;
     }

    return result;
}
