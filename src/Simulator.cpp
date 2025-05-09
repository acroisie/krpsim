#include "Simulator.hpp"
#include <algorithm>
#include <iostream>
#include <queue>
#include <set>

using namespace std;

Simulator::Simulator(const Config &configRef, int timeLimit)
    : config(configRef), timeLimit(timeLimit) {
    for (const auto &process : config.getProcesses()) {
        processMap[process.name] = &process;
    }
    buildProcessPriority();
}

void Simulator::buildProcessPriority() {
    set<string> goalSet(config.getOptimizeGoal().begin(),
                        config.getOptimizeGoal().end());

    for (const auto &process : config.getProcesses()) {
        for (const auto &[resource, _] : process.outputs) {
            if (goalSet.count(resource)) {
                processPriority[process.name] = 0;
                break;
            }
        }
    }
    for (int depth = 1; depth <= 2; ++depth) {
        for (const auto &process : config.getProcesses()) {
            if (processPriority.count(process.name)) continue;
            bool useful = false;
            for (const auto &[outputResource, _] : process.outputs) {
                for (const auto &otherProcess : config.getProcesses()) {
                    if (processPriority.count(otherProcess.name) &&
                        processPriority[otherProcess.name] == depth - 1 &&
                        otherProcess.inputs.count(outputResource)) {
                        useful = true;
                        break;
                    }
                }
                if (useful) break;
            }
            if (useful) processPriority[process.name] = depth;
        }
    }
}

const Process *Simulator::getProcessByName(const string &name) const {
    auto processIt = processMap.find(name);
    return processIt != processMap.end() ? processIt->second : nullptr;
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
    unordered_map<string, int> initialStocks;
    for (const auto &stock : config.getStocks()) {
        result.finalStocks[stock.name] = stock.quantity;
        initialStocks[stock.name] = stock.quantity;
    }

    size_t processSequenceIndex = 0;
    bool opportunistMode = false;
    int currentTime = 0;
    int processesExecuted = 0;

    priority_queue<RunningProcess, vector<RunningProcess>,
                   greater<RunningProcess>>
        runningProcessesQueue;

    while (currentTime <= timeLimit) {
        bool stocksUpdated = false;
        bool startedProcess = false;

        while (!runningProcessesQueue.empty() &&
               runningProcessesQueue.top().completionTime <= currentTime) {
            RunningProcess finishedProcess = runningProcessesQueue.top();
            runningProcessesQueue.pop();
            if (finishedProcess.processPtr) {
                for (const auto &[resource, quantity] :
                     finishedProcess.processPtr->outputs) {
                    result.finalStocks[resource] += quantity;
                }
                ++processesExecuted;
                stocksUpdated = true;
            }
        }
        if (currentTime >= timeLimit) {
            result.timeoutReached = true;
            break;
        }

        bool tryMore = true;
        while (tryMore) {
            tryMore = false;
            const Process *chosenProcess = nullptr;
            string chosenProcessName;
            bool fromSequence = false;

            if (!opportunistMode &&
                processSequenceIndex < processSequence.size()) {
                chosenProcessName = processSequence[processSequenceIndex];
                chosenProcess = getProcessByName(chosenProcessName);
                fromSequence = true;
            } else {
                opportunistMode = true;
                vector<const Process *> candidateProcesses;
                candidateProcesses.reserve(config.getProcesses().size());
                for (const auto &process : config.getProcesses())
                    candidateProcesses.push_back(&process);
                sort(candidateProcesses.begin(), candidateProcesses.end(),
                     [&](const Process *processA, const Process *processB) {
                         int priorityA = processPriority.count(processA->name)
                                             ? processPriority[processA->name]
                                             : 3;
                         int priorityB = processPriority.count(processB->name)
                                             ? processPriority[processB->name]
                                             : 3;
                         if (priorityA != priorityB)
                             return priorityA < priorityB;
                         return processA->nbCycle < processB->nbCycle;
                     });
                for (const Process *process : candidateProcesses) {
                    if (canStartProcess(process, result.finalStocks)) {
                        chosenProcess = process;
                        chosenProcessName = process->name;
                        break;
                    }
                }
            }

            if (!chosenProcess) {
                if (fromSequence) ++processSequenceIndex;
                continue;
            }
            if (!canStartProcess(chosenProcess, result.finalStocks) ||
                currentTime + chosenProcess->nbCycle > timeLimit) {
                if (fromSequence) ++processSequenceIndex;
                continue;
            }

            for (const auto &[resource, quantity] : chosenProcess->inputs) {
                result.finalStocks[resource] -= quantity;
            }
            runningProcessesQueue.push(
                {chosenProcess, currentTime + chosenProcess->nbCycle});
            result.executionLog.emplace_back(currentTime, chosenProcessName);
            startedProcess = true;
            tryMore = true;
            if (fromSequence) ++processSequenceIndex;
        }

        if (!startedProcess && !stocksUpdated) {
            if (runningProcessesQueue.empty()) break;
            currentTime =
                min(runningProcessesQueue.top().completionTime, timeLimit);
        } else {
            ++currentTime;
        }
    }

    int lastCycle = 0;
    for (const auto &[startCycle, processName] : result.executionLog) {
        const Process *processPtr = getProcessByName(processName);
        if (processPtr)
            lastCycle = max(lastCycle, startCycle + processPtr->nbCycle);
    }
    while (!runningProcessesQueue.empty() &&
           runningProcessesQueue.top().completionTime <= lastCycle) {
        const Process *processPtr = runningProcessesQueue.top().processPtr;
        runningProcessesQueue.pop();
        if (processPtr) {
            for (const auto &[resource, quantity] : processPtr->outputs)
                result.finalStocks[resource] += quantity;
        }
    }
    result.finalCycle = lastCycle;

    const auto &optimizationGoals = config.getOptimizeGoal();
    bool optimizeTime = find(optimizationGoals.begin(), optimizationGoals.end(),
                             "time") != optimizationGoals.end();

    result.fitness = 0.0;
    for (const string &goal : optimizationGoals) {
        if (goal == "time") continue;
        long delta = result.finalStocks[goal] - initialStocks[goal];
        result.fitness += static_cast<double>(delta);
    }
    if (optimizeTime) result.fitness /= (1.0 + result.finalCycle);

    result.fitness += processesExecuted * (optimizeTime ? 1e-3 : 1e-2);

    if (result.fitness == 0.0 && processesExecuted == 0) result.fitness = -1e9;

    return result;
}