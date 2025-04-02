#include "Simulator.hpp"
#include <iostream> // Pour debug éventuel
#include <limits>   // Pour std::numeric_limits
#include <queue>    // Pour std::priority_queue
#include <map>
#include <vector>
#include <string>
#include <set>      // Pour suivre les processus déjà tentés
#include <algorithm> // Pour std::min
#include <cmath>     // Pour std::isfinite

// Constructeur du Simulateur
Simulator::Simulator(const Config& config, int timeLimit)
    : systemConfig_(config), simulationTimeLimit_(timeLimit) {}

// Helper : Trouve un processus par son nom
const Process* Simulator::findProcessByName(const std::string& name) const {
    for (const auto& proc : systemConfig_.getProcesses()) {
        if (proc.name == name) {
            return &proc;
        }
    }
    return nullptr;
}

// Helper : Vérifie si un processus peut démarrer avec les stocks actuels
bool Simulator::canProcessStart(const Process* process, const std::map<std::string, int>& currentStocks) const {
    if (!process) {
        return false;
    }
    for (const auto& input : process->inputs) {
        auto it = currentStocks.find(input.first);
        if (it == currentStocks.end() || it->second < input.second) {
            return false;
        }
    }
    return true;
}

// Fonction CLÉ : Calcule le fitness et génère les logs d'exécution
double Simulator::calculateFitnessAndLogs(
    const std::vector<std::string>& sequenceToAttempt,
    std::vector<std::pair<int, std::string>>& executionLogs) {

    // --- 1. Initialisation de l'état de la simulation ---
    std::map<std::string, int> currentStocks;
    for (const auto& stock : systemConfig_.getStocks()) {
        currentStocks[stock.name] = stock.quantity;
    }
    executionLogs.clear();
    int currentTime = 0;
    int executedProcessCount = 0;
    size_t sequenceIndex = 0;

    std::priority_queue<RunningProcessInfo, std::vector<RunningProcessInfo>, std::greater<RunningProcessInfo>> runningProcessesQueue;
    bool opportunisticMode = false;

    // --- 2. Boucle de Simulation Principale ---
    while (currentTime < simulationTimeLimit_) {

        // --- A. Compléter les processus terminés ---
        bool resourcesUpdated = false;
        while (!runningProcessesQueue.empty() && runningProcessesQueue.top().completionTime <= currentTime) {
            RunningProcessInfo finishedProcess = runningProcessesQueue.top();
            runningProcessesQueue.pop();
            if (finishedProcess.processDetails) {
                 for (const auto& output : finishedProcess.processDetails->outputs) {
                    currentStocks[output.first] += output.second;
                    resourcesUpdated = true;
                 }
            }
        }

        // --- B. Essayer de lancer de nouveaux processus ---
        bool processStartedThisCycle = false;
        bool tryStartingMore = true;
        while (tryStartingMore) {
            tryStartingMore = false;
            const Process* processToTry = nullptr;
            std::string processNameToLog;
            bool attemptedFromSequence = false;

            // Choisir quel processus essayer:
            if (!opportunisticMode && sequenceIndex < sequenceToAttempt.size()) {
                // Mode guidé: essayer le prochain de la séquence
                processNameToLog = sequenceToAttempt[sequenceIndex];
                processToTry = findProcessByName(processNameToLog);
                attemptedFromSequence = true;
                // NE PAS incrémenter sequenceIndex ici, seulement si succès
            } else {
                 // Mode opportuniste
                 opportunisticMode = true;
                 const Process* bestOpportunisticChoice = nullptr;
                 for(const auto& potentialProcess : systemConfig_.getProcesses()) {
                     if (canProcessStart(&potentialProcess, currentStocks)) {
                         bestOpportunisticChoice = &potentialProcess;
                         processNameToLog = potentialProcess.name;
                         break;
                     }
                 }
                 processToTry = bestOpportunisticChoice;
                 if (!processToTry) break; // Rien à faire en opportuniste
            }

            // Si on a trouvé un processus à essayer
            if (processToTry) {
                if (canProcessStart(processToTry, currentStocks)) {
                    // Lancer le processus
                    for (const auto& input : processToTry->inputs) {
                        currentStocks[input.first] -= input.second;
                    }
                    runningProcessesQueue.push({processToTry, currentTime + processToTry->nbCycle});
                    executionLogs.push_back({currentTime, processNameToLog});
                    executedProcessCount++;
                    processStartedThisCycle = true;
                    tryStartingMore = true; // Essayer d'en lancer un autre

                    // Si on a réussi en mode guidé, on incrémente l'index de la séquence
                    if (attemptedFromSequence) {
                        sequenceIndex++;
                    }
                } else {
                    // Échec du lancement de CE processus
                    // Si on était en mode guidé, on passe au suivant dans la séquence lors de la prochaine itération externe
                     if (attemptedFromSequence) {
                         sequenceIndex++; // On passe au suivant de la séquence même si celui-ci a échoué
                     }
                     // Si on était en opportuniste et qu'on échoue (ne devrait pas arriver si canProcessStart est appelé avant), on break
                }
            } else {
                 // processToTry est nullptr (nom invalide ou fin de séquence)
                 if (attemptedFromSequence) {
                      sequenceIndex++; // Passer outre le nom invalide dans la séquence
                 }
                 // Si on est en fin de séquence, on passera en mode opportuniste au tour suivant
                 if (!opportunisticMode && sequenceIndex >= sequenceToAttempt.size()) {
                     opportunisticMode = true;
                     tryStartingMore = true; // Relancer la boucle interne pour tester l'opportuniste
                 }
            }
        } // Fin de la boucle interne while(tryStartingMore)

        // --- C. Avancer le Temps ---
        if (processStartedThisCycle || resourcesUpdated) {
            continue; // Rester au même temps pour voir si autre chose peut démarrer
        } else {
            // Rien n'a démarré, rien ne s'est terminé
            if (runningProcessesQueue.empty()) {
                // Plus rien en cours et rien ne peut démarrer -> Simulation terminée/bloquée
                break;
            } else {
                // Avancer au temps de complétion du prochain processus
                int nextCompletionTime = runningProcessesQueue.top().completionTime;
                if (nextCompletionTime >= simulationTimeLimit_) {
                    currentTime = simulationTimeLimit_;
                    break;
                }
                currentTime = nextCompletionTime;
            }
        }
    } // Fin de la boucle while (currentTime < simulationTimeLimit_)

    // --- D. Nettoyage final ---
     currentTime = std::min(currentTime, simulationTimeLimit_);

    // --- 3. Calcul Final du Fitness ---
    double fitness = 0.0;
    const auto& optimizationGoals = systemConfig_.getOptimizeGoal();

    // *** CORRECTION ICI ***
    // Retourner un grand négatif au lieu de lowest() si rien n'est exécuté
    if (executedProcessCount == 0) {
        // Utiliser une valeur négative importante mais évitable pour les calculs de shift
        // std::cerr << "Warning: No processes executed for a sequence." << std::endl; // Debug
        return -1.0e18;
    }

    // Calcul basé sur les objectifs
    bool timeIsGoal = false;
    int finalEventTime = 0; // Calculer le temps final pour l'objectif 'time'

    if (!executionLogs.empty()) {
        const auto& lastLog = executionLogs.back();
        const Process* lastProc = findProcessByName(lastLog.second);
        if(lastProc) {
            finalEventTime = lastLog.first + lastProc->nbCycle;
        }
    }
    finalEventTime = std::min(finalEventTime, simulationTimeLimit_);


    for (const std::string& goal : optimizationGoals) {
        if (goal == "time") {
            timeIsGoal = true;
            // Fitness inversement proportionnel au temps.
            fitness += 10000.0 / (1.0 + finalEventTime);
        } else {
            // Maximiser un stock spécifique
            if (currentStocks.count(goal)) {
                fitness += static_cast<double>(currentStocks.at(goal));
            }
        }
    }

    // Bonus/Malus additionnels:
    fitness += executedProcessCount * 0.01; // Bonus activité

    // Malus si bloqué avant la fin (et time n'est pas l'objectif)
    if (!timeIsGoal && currentTime < simulationTimeLimit_ && runningProcessesQueue.empty() && opportunisticMode) {
         fitness *= 0.5; // Pénalité blocage
    }

    // Vérifier si le fitness est valide avant de retourner
     if (!std::isfinite(fitness)) {
          std::cerr << "Warning: Calculated fitness is not finite (" << fitness << "). Resetting to large negative." << std::endl;
          return -1.0e17; // Retourner une autre grande valeur négative
     }


    return fitness;
}
