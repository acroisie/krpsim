#include "Simulator.hpp"
#include <algorithm>
#include <iostream>
#include <queue>
#include <set>

using namespace std;

using ProcessPriority = std::unordered_map<std::string, int>;

Simulator::Simulator(const Config &cfg, int limit)
    : config(cfg), timeLimit(limit) {
    for (const auto &p : config.getProcesses()) processMap[p.name] = &p;
    buildProcessPriority();
}

void Simulator::buildProcessPriority() {
    // Build the process priority table: 0 = direct goal, 1/2 = indirect
    std::set<std::string> goalNames(config.getOptimizeGoal().begin(),
                                    config.getOptimizeGoal().end());
    for (const auto &process : config.getProcesses())
        for (const auto &[outputName, _] : process.outputs)
            if (goalNames.count(outputName)) processPriority[process.name] = 0;
    for (int depth = 1; depth <= 2; ++depth)
        for (const auto &process : config.getProcesses()) {
            if (processPriority.count(process.name)) continue;
            for (const auto &[outputName, _] : process.outputs)
                for (const auto &otherProcess : config.getProcesses())
                    if (processPriority.count(otherProcess.name) &&
                        processPriority[otherProcess.name] == depth - 1 &&
                        otherProcess.inputs.count(outputName))
                        processPriority[process.name] = depth;
        }
}

const Process *
Simulator::getProcessByName(const std::string &processName) const {
    auto it = processMap.find(processName);
    return it != processMap.end() ? it->second : nullptr;
}

bool Simulator::canStartProcess(
    const Process *process,
    const std::map<std::string, int> &availableStocks) const {
    if (!process) return false;
    for (const auto &[resource, quantity] : process->inputs)
        if (availableStocks.count(resource) == 0 ||
            availableStocks.at(resource) < quantity)
            return false;
    return true;
}

// Pick the best process among candidates (priority then duration)
static const Process *
pickBestProcess(const std::vector<const Process *> &candidates,
                const ProcessPriority &priorityTable) {
    return *std::min_element(
        candidates.begin(), candidates.end(),
        [&](const Process *a, const Process *b) {
            int pa =
                priorityTable.count(a->name) ? priorityTable.at(a->name) : 3;
            int pb =
                priorityTable.count(b->name) ? priorityTable.at(b->name) : 3;
            return pa != pb ? pa < pb : a->nbCycle < b->nbCycle;
        });
}

Simulator::SimulationResult
Simulator::runSimulation(const std::vector<std::string> &processSequence) {
    SimulationResult result;
    std::unordered_map<std::string, int> initialStocks;
    for (const auto &stock : config.getStocks())
        result.finalStocks[stock.name] = initialStocks[stock.name] =
            stock.quantity;
    size_t sequenceIndex = 0;
    bool opportunistMode = false;
    int currentTime = 0, executedProcessCount = 0;
    std::priority_queue<RunningProcess, std::vector<RunningProcess>,
                        std::greater<RunningProcess>>
        runningProcesses;
    while (currentTime <= timeLimit) {
        bool stocksUpdated = false, startedProcess = false;
        // Complete finished processes
        while (!runningProcesses.empty() &&
               runningProcesses.top().completionTime <= currentTime) {
            auto finished = runningProcesses.top();
            runningProcesses.pop();
            if (finished.processPtr)
                for (const auto &[resource, quantity] :
                     finished.processPtr->outputs)
                    result.finalStocks[resource] += quantity;
            ++executedProcessCount;
            stocksUpdated = true;
        }
        if (currentTime >= timeLimit) {
            result.timeoutReached = true;
            break;
        }
        // Try to start as many processes as possible
        bool tryMore = true;
        while (tryMore) {
            tryMore = false;
            const Process *chosenProcess = nullptr;
            std::string chosenName;
            bool fromSequence = false;
            if (!opportunistMode && sequenceIndex < processSequence.size()) {
                chosenName = processSequence[sequenceIndex];
                chosenProcess = getProcessByName(chosenName);
                fromSequence = true;
            } else {
                opportunistMode = true;
                std::vector<const Process *> candidates;
                for (const auto &process : config.getProcesses())
                    if (canStartProcess(&process, result.finalStocks))
                        candidates.push_back(&process);
                if (!candidates.empty()) {
                    chosenProcess =
                        pickBestProcess(candidates, processPriority);
                    chosenName = chosenProcess->name;
                }
            }
            if (!chosenProcess) {
                if (fromSequence) ++sequenceIndex;
                continue;
            }
            if (!canStartProcess(chosenProcess, result.finalStocks) ||
                currentTime + chosenProcess->nbCycle > timeLimit) {
                if (fromSequence) ++sequenceIndex;
                continue;
            }
            for (const auto &[resource, quantity] : chosenProcess->inputs)
                result.finalStocks[resource] -= quantity;
            runningProcesses.push(
                {chosenProcess, currentTime + chosenProcess->nbCycle});
            result.executionLog.emplace_back(currentTime, chosenName);
            startedProcess = true;
            tryMore = true;
            if (fromSequence) ++sequenceIndex;
        }
        if (!startedProcess && !stocksUpdated) {
            if (runningProcesses.empty()) break;
            currentTime =
                std::min(runningProcesses.top().completionTime, timeLimit);
        } else {
            ++currentTime;
        }
    }
    // Finalize all processes that finish before the time limit
    while (!runningProcesses.empty() &&
           runningProcesses.top().completionTime <= timeLimit) {
        const Process *process = runningProcesses.top().processPtr;
        runningProcesses.pop();
        if (process)
            for (const auto &[resource, quantity] : process->outputs)
                result.finalStocks[resource] += quantity;
    }
    // Compute the actual duration
    for (const auto &[start, name] : result.executionLog) {
        const Process *process = getProcessByName(name);
        if (process)
            result.finalCycle =
                std::max(result.finalCycle, start + process->nbCycle);
    }
    result.finalCycle = std::min(result.finalCycle, timeLimit);
    // Compute the fitness score
    const auto &goals = config.getOptimizeGoal();
    bool optimizeTime =
        std::find(goals.begin(), goals.end(), "time") != goals.end();
    result.fitness = 0.0;
    for (const std::string &goal : goals)
        if (goal != "time")
            result.fitness += result.finalStocks[goal] - initialStocks[goal];
    if (optimizeTime) result.fitness /= (1.0 + result.finalCycle);
    // Tie-breaker: reward more executed processes
    result.fitness += executedProcessCount * (optimizeTime ? 1e-3 : 1e-2);
    if (result.fitness == 0.0 && executedProcessCount == 0)
        result.fitness = -1e9;
    return result;
}