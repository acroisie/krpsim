#include "Simulator.hpp"
#include <algorithm>
#include <iostream>
#include <queue>
#include <set>

using namespace std;

namespace {
    void applyProcessOutputs(const Process *process,
                             std::map<string, int> &stocks) {
        if (!process) return;
        for (const auto &[resource, quantity] : process->outputs)
            stocks[resource] += quantity;
    }
    int findLastCycle(const vector<pair<int, string>> &log,
                      const Simulator *sim) {
        int last = 0;
        for (const auto &[start, name] : log) {
            const Process *p = sim->getProcessByName(name);
            if (p) last = max(last, start + p->nbCycle);
        }
        return last;
    }
} // namespace

Simulator::Simulator(const Config &configRef, int timeLimit)
    : config(configRef), timeLimit(timeLimit) {
    for (const auto &process : config.getProcesses())
        processMap[process.name] = &process;
    buildProcessPriority();
}

void Simulator::buildProcessPriority() {
    set<string> goalSet(config.getOptimizeGoal().begin(),
                        config.getOptimizeGoal().end());
    for (const auto &process : config.getProcesses())
        for (const auto &[resource, _] : process.outputs)
            if (goalSet.count(resource)) {
                processPriority[process.name] = 0;
                break;
            }
    for (int depth = 1; depth <= 2; ++depth)
        for (const auto &process : config.getProcesses()) {
            if (processPriority.count(process.name)) continue;
            bool useful = false;
            for (const auto &[outputResource, _] : process.outputs)
                for (const auto &otherProcess : config.getProcesses())
                    if (processPriority.count(otherProcess.name) &&
                        processPriority[otherProcess.name] == depth - 1 &&
                        otherProcess.inputs.count(outputResource)) {
                        useful = true;
                        break;
                    }
            if (useful) processPriority[process.name] = depth;
        }
}

const Process *Simulator::getProcessByName(const string &name) const {
    auto it = processMap.find(name);
    return it != processMap.end() ? it->second : nullptr;
}

bool Simulator::canStartProcess(const Process *processPtr,
                                const map<string, int> &stocks) const {
    if (!processPtr) return false;
    for (const auto &[resource, quantity] : processPtr->inputs) {
        auto stockIt = stocks.find(resource);
        if (stockIt == stocks.end() || stockIt->second < quantity) return false;
    }
    return true;
}

Simulator::SimulationResult
Simulator::runSimulation(const vector<string> &processSequence) {
    SimulationResult result;
    for (const auto &stock : config.getStocks()) {
        result.finalStocks[stock.name] = stock.quantity;
    }
    map<string, int> initialStocks = result.finalStocks;
    size_t seqIdx = 0;
    bool opportunist = false;
    int currentTime = 0, processesExecuted = 0;
    priority_queue<RunningProcess, vector<RunningProcess>,
                   greater<RunningProcess>>
        runningQ;
    while (currentTime <= timeLimit) {
        bool stocksUpdated = false, startedProcess = false;
        while (!runningQ.empty() &&
               runningQ.top().completionTime <= currentTime) {
            applyProcessOutputs(runningQ.top().processPtr, result.finalStocks);
            runningQ.pop();
            ++processesExecuted;
            stocksUpdated = true;
        }
        if (currentTime >= timeLimit) {
            result.timeoutReached = true;
            break;
        }
        bool tryMore = true;
        while (tryMore) {
            tryMore = false;
            const Process *chosen = nullptr;
            string chosenName;
            bool fromSeq = false;
            if (!opportunist && seqIdx < processSequence.size()) {
                chosenName = processSequence[seqIdx];
                chosen = getProcessByName(chosenName);
                fromSeq = true;
            } else {
                opportunist = true;
                vector<const Process *> candidates;
                for (const auto &p : config.getProcesses())
                    candidates.push_back(&p);
                sort(candidates.begin(), candidates.end(),
                     [&](const Process *a, const Process *b) {
                         int pa = processPriority.count(a->name)
                                      ? processPriority[a->name]
                                      : 3;
                         int pb = processPriority.count(b->name)
                                      ? processPriority[b->name]
                                      : 3;
                         return pa != pb ? pa < pb : a->nbCycle < b->nbCycle;
                     });
                for (const Process *p : candidates)
                    if (canStartProcess(p, result.finalStocks)) {
                        chosen = p;
                        chosenName = p->name;
                        break;
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
            for (const auto &[resource, quantity] : chosen->inputs)
                result.finalStocks[resource] -= quantity;
            runningQ.push({chosen, currentTime + chosen->nbCycle});
            result.executionLog.emplace_back(currentTime, chosenName);
            startedProcess = true;
            tryMore = true;
            if (fromSeq) ++seqIdx;
        }
        if (!startedProcess && !stocksUpdated) {
            if (runningQ.empty()) break;
            currentTime = min(runningQ.top().completionTime, timeLimit);
        } else
            ++currentTime;
    }
    int lastCycle = findLastCycle(result.executionLog, this);
    while (!runningQ.empty() && runningQ.top().completionTime <= lastCycle) {
        applyProcessOutputs(runningQ.top().processPtr, result.finalStocks);
        runningQ.pop();
    }
    result.finalCycle = lastCycle;
    const auto &goals = config.getOptimizeGoal();
    bool optimiseTime = find(goals.begin(), goals.end(), "time") != goals.end();
    result.fitness = 0.0;
    for (const string &g : goals) {
        if (g == "time") continue;
        result.fitness +=
            static_cast<double>(result.finalStocks[g] - initialStocks[g]);
    }
    if (optimiseTime) result.fitness /= (1.0 + result.finalCycle);
    result.fitness += processesExecuted * (optimiseTime ? 1e-3 : 1e-2);
    if (result.fitness == 0.0 && processesExecuted == 0) result.fitness = -1e9;
    return result;
}