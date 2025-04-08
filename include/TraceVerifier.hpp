#pragma once
#include "Config.hpp"
#include "Process.hpp"
#include <queue>
#include <map>
#include <set>
#include <string>

/**
 *  TraceVerifier
 *  -------------
 *  Vérifie qu’un fichier-trace “<cycle>:<process>” est cohérent avec la
 *  configuration krpsim : disponibilité des stocks au démarrage, durées
 *  d’exécution, etc.  Affiche l’erreur précise ou le stock final.
 */
class TraceVerifier {
  public:
    TraceVerifier(const Config &cfg, int timeLimit = 1'000'000);
    /** Renvoie true si la trace est valide, false sinon (affiche les détails). */
    bool verifyFile(const std::string &traceFilename);

  private:
    const Config &config;
    int timeLimit;
    std::map<std::string, int> stocks;

    struct RunningProcess {
        const Process *proc;
        int completionTime;
        bool operator>(const RunningProcess &o) const {
            return completionTime > o.completionTime;
        }
    };
    std::priority_queue<RunningProcess,
                        std::vector<RunningProcess>,
                        std::greater<RunningProcess>>
        runningQ;

    bool canStart(const Process *p) const;
    void finishUntil(int cycle);
    bool start(int cycle, const Process *p, std::string &err);
    void addOutputs(const Process *p);
    void printStocks(int finalCycle) const;
};
