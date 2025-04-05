#pragma once
#include "Config.hpp"
#include "Process.hpp"
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

/**
 *  Simulator
 *  ---------
 *  Exécute une séquence de noms de process et – en mode opportuniste –
 *  complète avec d'autres process classés par priorité afin de maximiser
 *  les objectifs déclarés dans la ligne `optimize:` du fichier krpsim.
 */
class Simulator {
  public:
    struct RunningProcess {
        const Process *processPtr;
        int completionTime;
        bool operator>(const RunningProcess &other) const {
            return completionTime > other.completionTime;
        }
    };

    struct SimulationResult {
        double fitness = std::numeric_limits<double>::lowest();
        std::vector<std::pair<int, std::string>> executionLog;
        std::map<std::string, int> finalStocks;
        int finalCycle = 0;
        bool timeoutReached = false;
    };

    Simulator(const Config &config, int timeLimit);

    SimulationResult runSimulation(const std::vector<std::string> &processSequence);

    /** Accès lecture à la table de priorité (utile à GeneticAlgorithm). */
    const std::unordered_map<std::string, int> &getProcessPriority() const {
        return processPriority;
    }

  private:
    const Config &config;
    int timeLimit;

    // Cache nom ➜ pointeur Process pour accès O(1)
    std::unordered_map<std::string, const Process *> processMap;
    // Priorité calculée une fois : 0 = produit direct d'un objectif, 1..n sinon
    std::unordered_map<std::string, int> processPriority;

    // helpers
    void buildProcessPriority();
    bool canStartProcess(const Process *process, const std::map<std::string, int> &stocks) const;
    const Process *getProcessByName(const std::string &name) const;
};
