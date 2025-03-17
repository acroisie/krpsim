#pragma once

#include "Config.hpp"
#include <map>
#include <vector>
#include <string>
#include <queue>
#include <memory>
#include <functional>

class Simulator {
public:
    struct Event {
        int time;
        std::string processName;
        bool isStart;  // true = start, false = finish

        bool operator>(const Event& other) const {
            return time > other.time;
        }
    };

    struct Result {
        std::vector<std::pair<int, std::string>> executionLog;
        std::map<std::string, int> finalStocks;
        int finalTime;
        double score;
    };

    Simulator(const Config& config);

    // Simule l'exécution avec une priorité donnée pour chaque processus
    Result simulate(const std::vector<std::string>& processSequence, 
                    int timeLimit, 
                    bool generateLog = true);

    // Calcule un score en fonction des objectifs d'optimisation
    double calculateScore(const std::map<std::string, int>& stocks, int time) const;

private:
    const Config& config_;
    std::priority_queue<Event, std::vector<Event>, std::greater<Event>> eventQueue_;
    std::map<std::string, int> stocks_;
    std::vector<std::pair<int, std::string>> executionLog_;
    int currentTime_;

    bool canRunProcess(const Process& process) const;
    void startProcess(const Process& process, int startTime);
    void consumeResources(const Process& process);
    void produceResources(const Process& process);
    const Process* getProcessByName(const std::string& name) const;
};