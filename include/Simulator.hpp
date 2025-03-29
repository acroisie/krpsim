#pragma once
#include "Config.hpp"
#include "Process.hpp"
#include <map>
#include <string>
#include <vector>
#include <utility>

class Simulator {
public:
    Simulator(const Config& config, int delayLimit);
    
    // Simule l'exécution d'une séquence de processus
    bool simulateSequence(const std::vector<std::string>& sequence,
                          std::map<std::string, int>& finalStocks,
                          int& executedCount);
    
    // Calcule le fitness pour une séquence de processus
    double calculateFitness(const std::vector<std::string>& processSequence,
                          std::vector<std::pair<int, std::string>>& executionLogs);
    
    // Obtient tous les processus exécutables avec les stocks courants
    std::vector<const Process*> getRunnableProcesses(const std::map<std::string, int>& stocks);
    
    // Exécute un processus, modifiant les niveaux de stock
    bool executeProcess(const Process* process, std::map<std::string, int>& stocks);

private:
    const Config& config_;
    int delayLimit_;
};