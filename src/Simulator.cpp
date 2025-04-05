#include "Simulator.hpp"
#include <algorithm>
#include <iostream>
#include <queue>
#include <set>

using namespace std;

Simulator::Simulator(const Config &cfg, int limit) : config(cfg), timeLimit(limit) {
    for (const auto &p : config.getProcesses()) {
        processMap[p.name] = &p;
    }
    buildProcessPriority();
}

void Simulator::buildProcessPriority() {
    // 0 ➜ produit directement un objectif ; profondeur limitée à 2 suffisant ici.
    set<string> goalSet(config.getOptimizeGoal().begin(), config.getOptimizeGoal().end());

    for (const auto &p : config.getProcesses()) {
        for (const auto &[res, _] : p.outputs) {
            if (goalSet.count(res)) {
                processPriority[p.name] = 0;
                break;
            }
        }
    }
    for (int depth = 1; depth <= 2; ++depth) {
        for (const auto &p : config.getProcesses()) {
            if (processPriority.count(p.name)) continue;
            bool useful = false;
            for (const auto &[outRes, _] : p.outputs) {
                for (const auto &q : config.getProcesses()) {
                    if (processPriority.count(q.name) && processPriority[q.name] == depth - 1 &&
                        q.inputs.count(outRes)) {
                        useful = true;
                        break;
                    }
                }
                if (useful) break;
            }
            if (useful) processPriority[p.name] = depth;
        }
    }
}

const Process *Simulator::getProcessByName(const string &name) const {
    auto it = processMap.find(name);
    return it != processMap.end() ? it->second : nullptr;
}

bool Simulator::canStartProcess(const Process *proc, const map<string, int> &stocks) const {
    if (!proc) return false;
    for (const auto &[res, qty] : proc->inputs) {
        auto it = stocks.find(res);
        if (it == stocks.end() || it->second < qty) return false;
    }
    return true;
}

Simulator::SimulationResult Simulator::runSimulation(const vector<string> &sequence) {
    SimulationResult result;
    unordered_map<string, int> initialStocks;
    for (const auto &s : config.getStocks()) {
        result.finalStocks[s.name] = s.quantity;
        initialStocks[s.name]     = s.quantity;
    }

    size_t seqIdx           = 0;
    bool   opportunist      = false;
    int    currentTime      = 0;
    int    processesExecuted = 0;

    priority_queue<RunningProcess, vector<RunningProcess>, greater<RunningProcess>> runningQ;

    while (currentTime <= timeLimit) {
        bool stocksUpdated  = false;
        bool startedProcess = false;

        // Terminer les process achevés
        while (!runningQ.empty() && runningQ.top().completionTime <= currentTime) {
            RunningProcess done = runningQ.top();
            runningQ.pop();
            if (done.processPtr) {
                for (const auto &[res, qty] : done.processPtr->outputs) {
                    result.finalStocks[res] += qty;
                }
                ++processesExecuted;
                stocksUpdated = true;
            }
        }
        if (currentTime >= timeLimit) {
            result.timeoutReached = true;
            break;
        }

        // Démarrer de nouveaux process tant qu'on peut
        bool tryMore = true;
        while (tryMore) {
            tryMore = false;
            const Process *chosen = nullptr;
            string chosenName;
            bool fromSeq = false;

            if (!opportunist && seqIdx < sequence.size()) {
                chosenName = sequence[seqIdx];
                chosen     = getProcessByName(chosenName);
                fromSeq    = true;
            } else {
                opportunist = true;
                vector<const Process *> cand;
                cand.reserve(config.getProcesses().size());
                for (const auto &p : config.getProcesses()) cand.push_back(&p);
                sort(cand.begin(), cand.end(), [&](const Process *a, const Process *b) {
                    int pa = processPriority.count(a->name) ? processPriority[a->name] : 3;
                    int pb = processPriority.count(b->name) ? processPriority[b->name] : 3;
                    if (pa != pb) return pa < pb;
                    return a->nbCycle < b->nbCycle;
                });
                for (const Process *p : cand) {
                    if (canStartProcess(p, result.finalStocks)) {
                        chosen     = p;
                        chosenName = p->name;
                        break;
                    }
                }
            }

            if (!chosen) {
                if (fromSeq) ++seqIdx;
                continue;
            }
            if (!canStartProcess(chosen, result.finalStocks) ||
                currentTime + chosen->nbCycle > timeLimit) {
                if (fromSeq) ++seqIdx;
                continue;
            }

            for (const auto &[res, qty] : chosen->inputs) {
                result.finalStocks[res] -= qty;
            }
            runningQ.push({chosen, currentTime + chosen->nbCycle});
            result.executionLog.emplace_back(currentTime, chosenName);
            startedProcess = true;
            tryMore        = true;
            if (fromSeq) ++seqIdx;
        }

        if (!startedProcess && !stocksUpdated) {
            if (runningQ.empty()) break;
            currentTime = min(runningQ.top().completionTime, timeLimit);
        } else {
            ++currentTime;
        }
    }

    // Finaliser ce qui termine avant timeLimit
    while (!runningQ.empty() && runningQ.top().completionTime <= timeLimit) {
        const Process *p = runningQ.top().processPtr;
        runningQ.pop();
        if (p) {
            for (const auto &[res, qty] : p->outputs) result.finalStocks[res] += qty;
        }
    }

    // Durée réelle
    for (const auto &[start, name] : result.executionLog) {
        const Process *p = getProcessByName(name);
        if (p) result.finalCycle = max(result.finalCycle, start + p->nbCycle);
    }
    result.finalCycle = min(result.finalCycle, timeLimit);

    // Fitness
    const auto &goals = config.getOptimizeGoal();
    bool optimiseTime = find(goals.begin(), goals.end(), "time") != goals.end();

    result.fitness = 0.0;
    for (const string &g : goals) {
        if (g == "time") continue;
        long delta = result.finalStocks[g] - initialStocks[g];
        result.fitness += static_cast<double>(delta);
    }
    if (optimiseTime) result.fitness /= (1.0 + result.finalCycle);

    // léger tie‑breaker
    result.fitness += processesExecuted * (optimiseTime ? 1e-3 : 1e-2);

    if (result.fitness == 0.0 && processesExecuted == 0) result.fitness = -1e9;

    return result;
}