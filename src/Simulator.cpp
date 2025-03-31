#include "Simulator.hpp"
#include <algorithm>
#include <iostream>
#include <limits>
#include <queue>

using namespace std;

Simulator::Simulator(const Config &config, int delayLimit) : config_(config), delayLimit_(delayLimit) {}

bool Simulator::simulateSequence(const vector<string> &sequence, map<string, int> &finalStocks, int &executedCount) {
    // Reset state
    map<string, int> stocks;
    for (const auto &stock : config_.getStocks()) {
        stocks[stock.name] = stock.quantity;
    }
    int cycle = 0;
    executedCount = 0;

    struct RunningProcess {
        const Process *process;
        int completionCycle;
    };
    vector<RunningProcess> runningProcesses;

    size_t idx = 0;
    while (cycle < delayLimit_ && idx < sequence.size()) {
        // Complete any finished processes
        auto it = runningProcesses.begin();
        while (it != runningProcesses.end()) {
            if (it->completionCycle <= cycle) {
                for (const auto &out : it->process->outputs) {
                    stocks[out.first] += out.second;
                }
                it = runningProcesses.erase(it);
            } else {
                ++it;
            }
        }

        // Try to execute the next process
        const string &procName = sequence[idx];
        const Process *proc = nullptr;

        // Find the process
        for (const auto &p : config_.getProcesses()) {
            if (p.name == procName) {
                proc = &p;
                break;
            }
        }

        if (!proc) {
            // Unknown process, skip
            idx++;
            continue;
        }

        // Check if we can execute it
        bool canExecute = true;
        for (const auto &input : proc->inputs) {
            if (stocks[input.first] < input.second) {
                canExecute = false;
                break;
            }
        }

        if (canExecute) {
            // Execute the process
            for (const auto &input : proc->inputs) {
                stocks[input.first] -= input.second;
            }
            runningProcesses.push_back({proc, cycle + proc->nbCycle});
            executedCount++;
            idx++;
        } else {
            // Can't execute, advance time if processes are running
            if (!runningProcesses.empty()) {
                int nextCompletion = delayLimit_;
                for (const auto &rp : runningProcesses) {
                    nextCompletion = min(nextCompletion, rp.completionCycle);
                }
                cycle = nextCompletion;
            } else {
                // No processes running, skip this one
                idx++;
            }
        }
    }

    // Complete any remaining processes
    while (!runningProcesses.empty() && cycle < delayLimit_) {
        int nextCompletion = delayLimit_;
        for (const auto &rp : runningProcesses) {
            nextCompletion = min(nextCompletion, rp.completionCycle);
        }
        cycle = nextCompletion;

        auto it = runningProcesses.begin();
        while (it != runningProcesses.end()) {
            if (it->completionCycle <= cycle) {
                for (const auto &out : it->process->outputs) {
                    stocks[out.first] += out.second;
                }
                it = runningProcesses.erase(it);
            } else {
                ++it;
            }
        }
    }

    finalStocks = stocks;
    return executedCount > 0;
}

// Calculer le fitness de manière générique sans hardcoder pour des problèmes spécifiques
double Simulator::calculateFitness(const vector<string> &processSequence, vector<pair<int, string>> &executionLogs) {
    // Reset state for simulation
    map<string, int> stocks;
    for (const auto &stock : config_.getStocks()) {
        stocks[stock.name] = stock.quantity;
    }
    executionLogs.clear();
    int currentCycle = 0;

    struct RunningProcess {
        const Process *process;
        int completionCycle;
    };
    vector<RunningProcess> runningProcesses;

    // Count executed processes
    int executed = 0;
    size_t index = 0;

    while (currentCycle < delayLimit_ && index < processSequence.size()) {
        // Complete finished processes
        auto it = runningProcesses.begin();
        while (it != runningProcesses.end()) {
            if (it->completionCycle <= currentCycle) {
                for (const auto &out : it->process->outputs) {
                    stocks[out.first] += out.second;
                }
                it = runningProcesses.erase(it);
            } else {
                ++it;
            }
        }

        // Try to execute the current process
        const string &procName = processSequence[index];
        const Process *chosen = nullptr;

        // Find the process
        for (const auto &proc : config_.getProcesses()) {
            if (proc.name == procName) {
                chosen = &proc;
                break;
            }
        }

        if (!chosen) {
            // Unknown process name, skip
            index++;
            continue;
        }

        // Check if we can execute it
        bool canExecute = true;
        for (const auto &input : chosen->inputs) {
            if (stocks.count(input.first) == 0 || stocks.at(input.first) < input.second) {
                canExecute = false;
                break;
            }
        }

        if (canExecute) {
            // Consume inputs
            for (const auto &input : chosen->inputs) {
                stocks[input.first] -= input.second;
            }

            runningProcesses.push_back({chosen, currentCycle + chosen->nbCycle});
            executionLogs.push_back({currentCycle, chosen->name});
            executed++;
            index++;
        } else {
            // Can't execute this process now
            if (!runningProcesses.empty()) {
                // Fast-forward to next process completion
                int nextCompletion = delayLimit_;
                for (auto &rp : runningProcesses) {
                    nextCompletion = min(nextCompletion, rp.completionCycle);
                }
                currentCycle = nextCompletion;
            } else {
                // No running processes, skip this one
                index++;
            }
        }
    }

    // Complete any remaining processes
    while (!runningProcesses.empty() && currentCycle < delayLimit_) {
        int nextCompletion = delayLimit_;
        for (auto &rp : runningProcesses) {
            nextCompletion = min(nextCompletion, rp.completionCycle);
        }
        currentCycle = nextCompletion;

        auto it = runningProcesses.begin();
        while (it != runningProcesses.end()) {
            if (it->completionCycle <= currentCycle) {
                for (const auto &out : it->process->outputs) {
                    stocks[out.first] += out.second;
                }
                it = runningProcesses.erase(it);
            } else {
                ++it;
            }
        }
    }

    // If nothing executed, return a very low fitness
    if (executed == 0) {
        return -1000.0;
    }

    // Calculate fitness based on optimization goals
    double fitness = 0.0;
    const auto &goals = config_.getOptimizeGoal();

    if (!goals.empty()) {
        // Gérer chaque objectif d'optimisation
        for (const auto &goal : goals) {
            if (goal == "time") {
                // Minimiser le temps - plus petit est meilleur
                // Normaliser pour que la valeur reste positive: 1000 - currentCycle
                fitness += max(0.0, 1000.0 - currentCycle);
            } else {
                // Maximiser la quantité d'un stock
                if (stocks.count(goal) > 0) {
                    fitness += stocks.at(goal);
                }
            }
        }

        // Récompense pour le nombre de processus exécutés et l'utilisation des ressources
        fitness += executed * 0.1; // Petit bonus pour l'exécution de nombreux processus
    }

    return fitness;
}
vector<const Process *> Simulator::getRunnableProcesses(const map<string, int> &stocks) {
    vector<const Process *> result;
    for (const auto &proc : config_.getProcesses()) {
        bool canRun = true;
        for (const auto &in : proc.inputs) {
            if (stocks.count(in.first) == 0 || stocks.at(in.first) < in.second) {
                canRun = false;
                break;
            }
        }
        if (canRun) {
            result.push_back(&proc);
        }
    }
    return result;
}

bool Simulator::executeProcess(const Process *process, map<string, int> &stocks) {
    if (!process) return false;

    // Verify we have enough resources
    for (const auto &in : process->inputs) {
        if (stocks.count(in.first) == 0 || stocks.at(in.first) < in.second) {
            return false;
        }
    }

    // Consume inputs
    for (const auto &in : process->inputs) {
        stocks[in.first] -= in.second;
    }

    // Produce outputs
    for (const auto &out : process->outputs) {
        stocks[out.first] += out.second;
    }

    return true;
}