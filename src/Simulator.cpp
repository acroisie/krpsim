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
          std::cerr << "DEBUG: Warning: Simulator initialized with no processes defined." << std::endl;
     }
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

    // std::cerr << "\nDEBUG: Starting Simulation..." << std::endl; // DEBUG

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

    // --- Boucle de Simulation (Logique v11 - sequenceIndex avance seulement si succès ET depuis séquence) ---
    while (currentTime <= simulationTimeLimit_) {
        // std::cerr << "DEBUG: Loop Start. currentTime = " << currentTime << std::endl; // DEBUG
        bool resourcesUpdated = false;
        while (!runningProcessesQueue.empty() && runningProcessesQueue.top().completionTime <= currentTime) {
            RunningProcessInfo finishedProcess = runningProcessesQueue.top();
            runningProcessesQueue.pop();
            if (finishedProcess.completionTime <= simulationTimeLimit_ && finishedProcess.processDetails) {
                 // std::cerr << "DEBUG: t=" << currentTime << ": Finishing process '" << finishedProcess.processDetails->name << "'" << std::endl; // DEBUG
                 for (const auto& output : finishedProcess.processDetails->outputs) {
                    result.finalStocks[output.first] += output.second;
                 }
                 resourcesUpdated = true;
            }
        }
        if (currentTime >= simulationTimeLimit_) { result.reachedTimeLimit = true; break; }

        bool processStartedThisCycle = false;
        if (currentTime < simulationTimeLimit_) {
            bool tryStartingMore = true;
            while (tryStartingMore) {
                tryStartingMore = false;
                const Process* processToTry = nullptr;
                std::string processNameToLog;
                bool attemptedFromSequence = false;
                if (!opportunisticMode && sequenceIndex < sequenceToAttempt.size()) {
                    processNameToLog = sequenceToAttempt[sequenceIndex];
                    processToTry = findProcessByName(processNameToLog);
                    attemptedFromSequence = true;
                } else {
                    opportunisticMode = true;
                    const Process* bestOpportunisticChoice = nullptr;
                    for(const auto& pair : processLookupMap_) {
                        const Process* potentialProcess = pair.second;
                        if (canProcessStart(potentialProcess, result.finalStocks)) {
                            bestOpportunisticChoice = potentialProcess;
                            processNameToLog = potentialProcess->name;
                            break;
                        }
                    }
                    processToTry = bestOpportunisticChoice;
                    if (!processToTry) break;
                }
                if (processToTry) {
                    if (canProcessStart(processToTry, result.finalStocks) && (currentTime + processToTry->nbCycle <= simulationTimeLimit_)) {
                        // std::cerr << "DEBUG: t=" << currentTime << ": Starting process '" << processNameToLog << "'" << std::endl; // DEBUG
                        for (const auto& input : processToTry->inputs) { result.finalStocks[input.first] -= input.second; }
                        runningProcessesQueue.push({processToTry, currentTime + processToTry->nbCycle});
                        result.logs.push_back({currentTime, processNameToLog});
                        executedProcessCount++;
                        processStartedThisCycle = true;
                        tryStartingMore = true;
                        if (attemptedFromSequence) { sequenceIndex++; } // Avance SEULEMENT si succès ET depuis séquence
                    } else {
                        // Échec
                        // Ne PAS avancer sequenceIndex si échec depuis séquence
                         if(!attemptedFromSequence) tryStartingMore = false;
                    }
                } else {
                    if (attemptedFromSequence) { sequenceIndex++; } // Avance si nom invalide
                    if (!opportunisticMode && sequenceIndex >= sequenceToAttempt.size()) { opportunisticMode = true; tryStartingMore = true; }
                }
            }
        }

        if (processStartedThisCycle || resourcesUpdated) {
             if (currentTime >= simulationTimeLimit_) { result.reachedTimeLimit = true; break; }
             continue;
        } else {
            if (runningProcessesQueue.empty()) {
                 // std::cerr << "DEBUG: t=" << currentTime << ": Queue empty and no action. Breaking loop." << std::endl; // DEBUG
                 break;
            } else {
                int nextCompletionTime = runningProcessesQueue.top().completionTime;
                if (nextCompletionTime > simulationTimeLimit_) { currentTime = simulationTimeLimit_; result.reachedTimeLimit = true; }
                else {
                    if (nextCompletionTime <= currentTime) { currentTime++; } // Sécurité anti-blocage minimale
                    else { currentTime = nextCompletionTime; }
                }
            }
        }
         if (currentTime > simulationTimeLimit_) { currentTime = simulationTimeLimit_; result.reachedTimeLimit = true; break; }
    } // Fin while

    // --- Dernière passe de complétion ---
     while (!runningProcessesQueue.empty() && runningProcessesQueue.top().completionTime <= simulationTimeLimit_) {
         RunningProcessInfo finishedProcess = runningProcessesQueue.top();
         runningProcessesQueue.pop();
         if (finishedProcess.processDetails) {
              // std::cerr << "DEBUG: t=" << currentTime << ": Finishing process '" << finishedProcess.processDetails->name << "' in final pass." << std::endl; // DEBUG
              for (const auto& output : finishedProcess.processDetails->outputs) {
                 result.finalStocks[output.first] += output.second;
              }
         }
     }

    // --- Calcul Final du Fitness (SIMPLIFIÉ) ---
    result.fitness = 0.0;
    const auto& optimizationGoals = systemConfig_.getOptimizeGoal();

    if (executedProcessCount == 0 && result.logs.empty()) {
        result.fitness = -1.0e18;
        result.finalCycle = currentTime;
        return result;
    }

    bool timeIsGoal = false;
    if (!result.logs.empty()) {
        int max_end_time = 0;
        for(const auto& log : result.logs) {
            const Process* proc = findProcessByName(log.second);
            if (proc) { max_end_time = std::max(max_end_time, log.first + proc->nbCycle); }
        }
        result.finalCycle = std::min(max_end_time, simulationTimeLimit_);
    } else {
         result.finalCycle = 0;
    }

    // *** FITNESS SIMPLIFIÉE ***
    for (const std::string& goal : optimizationGoals) {
        if (goal == "time") {
            timeIsGoal = true;
            result.fitness += 10000.0 / (1.0 + result.finalCycle);
        } else {
            // Donner un poids UNIQUEMENT à l'objectif final demandé
            double stockValue = 0;
            if (result.finalStocks.count(goal)) {
                stockValue = static_cast<double>(result.finalStocks.at(goal));
            }
            // Mettre un poids élevé pour l'objectif final, 0 pour le reste
            double weight = 1000.0; // Poids générique élevé pour l'objectif
            result.fitness += stockValue * weight;
        }
    }

    // Garder un bonus d'activité TRÈS faible pour départager les ex æquo
    result.fitness += executedProcessCount * 1e-9;

    // Pénalité de blocage (inchangée)
    if (!timeIsGoal && currentTime <= simulationTimeLimit_ && runningProcessesQueue.empty() && opportunisticMode) {
         result.fitness *= 0.5;
    }
     if (!std::isfinite(result.fitness)) {
          std::cerr << "Warning: Calculated fitness is not finite (" << result.fitness << "). Resetting to large negative." << std::endl;
          result.fitness = -1.0e17;
     }

    // std::cerr << "DEBUG: Final Fitness = " << result.fitness << std::endl; // DEBUG
    return result;
}
