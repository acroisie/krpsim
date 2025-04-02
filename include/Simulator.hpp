#pragma once

#include "Config.hpp"
#include "Process.hpp"
#include <map>
#include <string>
#include <vector>
#include <utility>
#include <queue>
#include <unordered_map>
#include <limits> // Pour std::numeric_limits

class Simulator {
  public:
    // Structure pour suivre un processus en cours (publique pour info si besoin externe)
    struct RunningProcessInfo {
        const Process* processDetails;
        int completionTime;
        bool operator>(const RunningProcessInfo& other) const {
            return completionTime > other.completionTime;
        }
    };

    // *** NOUVELLE STRUCTURE pour retourner les résultats complets ***
    struct SimulationResult {
        double fitness = std::numeric_limits<double>::lowest(); // Score calculé
        std::vector<std::pair<int, std::string>> logs;         // Log d'exécution <cycle, nom>
        std::map<std::string, int> finalStocks;                // État des stocks à la fin
        int finalCycle = 0;                                    // Cycle de fin de la simulation/dernier événement
        bool reachedTimeLimit = false;                         // Indique si la limite de temps a été atteinte
    };


  public:
    // Constructeur
    Simulator(const Config& config, int timeLimit);

    // *** MODIFIÉ : La fonction clé retourne maintenant SimulationResult ***
    SimulationResult calculateFitnessAndLogs(
        const std::vector<std::string>& sequenceToAttempt);

  private:
    const Config& systemConfig_;
    int simulationTimeLimit_;
    std::unordered_map<std::string, const Process*> processLookupMap_;

    // Helpers (inchangés)
    bool canProcessStart(const Process* process, const std::map<std::string, int>& currentStocks) const;
    const Process* findProcessByName(const std::string& name) const;
};
