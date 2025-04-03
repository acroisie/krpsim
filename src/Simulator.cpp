#include "Simulator.hpp"
#include <algorithm>
#include <iostream>
#include <set>

using namespace std;

Simulator::Simulator(const Config &config, int timeLimit)
    : config(config), timeLimit(timeLimit) {

    for (const auto &process : config.getProcesses()) {
        processMap[process.name] = &process;
    }
}

const Process *Simulator::getProcessByName(const std::string &name) const {
    auto it = processMap.find(name);
    if (it != processMap.end()) {
        return it->second;
    }
    return nullptr;
}

bool Simulator::canStartProcess(const Process *process,
                                const map<string, int> &stocks) const {
    if (!process) {
        return false;
    }

    for (const auto &[resource, quantity] : process->inputs) {
        auto stockIt = stocks.find(resource);
        if (stockIt == stocks.end() || stockIt->second < quantity) {
            return false;
        }
    }

    return true;
}

Simulator::SimulationResult
Simulator::runSimulation(const vector<string> &processSequence) {
    SimulationResult result;

    // Initialize stocks from config
    for (const auto &stock : config.getStocks()) {
        result.finalStocks[stock.name] = stock.quantity;
    }

    int currentTime = 0;
    int processesExecuted = 0;
    size_t sequenceIndex = 0;

    // Priority queue for running processes, ordered by completion time
    priority_queue<RunningProcess, vector<RunningProcess>,
                   greater<RunningProcess>>
        runningProcesses;
    bool opportunisticMode = false;
    
    // Keep track of resources that are currently in use by running processes
    map<string, int> resourcesInUse;

    // Map to track which process uses which resources (for better tracking)
    map<string, set<string>> processResourceMap;
    
    // Pre-compute which resources are reusable (same resource is both input and output)
    map<string, set<string>> reusableResources;
    for (const auto& process : config.getProcesses()) {
        set<string> reusable;
        for (const auto& [resource, quantity] : process.inputs) {
            auto outputIt = process.outputs.find(resource);
            if (outputIt != process.outputs.end() && outputIt->second >= quantity) {
                reusable.insert(resource);
            }
        }
        if (!reusable.empty()) {
            reusableResources[process.name] = reusable;
        }
    }

    // Main simulation loop
    while (currentTime <= timeLimit) {
        bool stocksUpdated = false;

        // Complete processes that have finished
        while (!runningProcesses.empty() &&
               runningProcesses.top().completionTime <= currentTime) {
            RunningProcess finishedProcess = runningProcesses.top();
            runningProcesses.pop();

            if (finishedProcess.completionTime <= timeLimit &&
                finishedProcess.processPtr) {
                for (const auto &[resource, quantity] :
                     finishedProcess.processPtr->outputs) {
                    result.finalStocks[resource] += quantity;
                }
                
                // Release resources that were being used by this process
                const string& processName = finishedProcess.processPtr->name;
                if (reusableResources.find(processName) != reusableResources.end()) {
                    for (const auto& resource : reusableResources[processName]) {
                        auto inputIt = finishedProcess.processPtr->inputs.find(resource);
                        if (inputIt != finishedProcess.processPtr->inputs.end()) {
                            resourcesInUse[resource] -= inputIt->second;
                            if (resourcesInUse[resource] <= 0) {
                                resourcesInUse.erase(resource);
                            }
                        }
                    }
                }
                
                stocksUpdated = true;
                processesExecuted++;
            }
        }

        if (currentTime >= timeLimit) {
            result.timeoutReached = true;
            break;
        }

        bool processStarted = false;
        
        // Try to start processes in a loop to allow multiple processes to start in the same cycle
        bool tryMoreProcesses = true;
        while (tryMoreProcesses && currentTime < timeLimit) {
            tryMoreProcesses = false;
            const Process *processToRun = nullptr;
            string processName;
            bool attemptedFromSequence = false;

            // Try to get a process from the sequence or opportunistically
            if (!opportunisticMode && sequenceIndex < processSequence.size()) {
                processName = processSequence[sequenceIndex];
                processToRun = getProcessByName(processName);
                attemptedFromSequence = true;
            } else {
                opportunisticMode = true;
                // Try all processes in order of config (deterministic)
                for (const auto &process : config.getProcesses()) {
                    if (canStartProcess(&process, result.finalStocks)) {
                        bool canUseResources = true;
                        // Check if reusable resources are available
                        if (reusableResources.find(process.name) != reusableResources.end()) {
                            for (const auto& resource : reusableResources[process.name]) {
                                auto inputIt = process.inputs.find(resource);
                                if (inputIt != process.inputs.end()) {
                                    auto inUseIt = resourcesInUse.find(resource);
                                    int inUseQty = (inUseIt != resourcesInUse.end()) ? inUseIt->second : 0;
                                    int availableQty = result.finalStocks[resource] - inUseQty;
                                    if (availableQty < inputIt->second) {
                                        canUseResources = false;
                                        break;
                                    }
                                }
                            }
                        }
                        
                        if (canUseResources) {
                            processToRun = &process;
                            processName = process.name;
                            break;
                        }
                    }
                }
            }

            if (processToRun) {
                // Check if resources are available and not in use
                bool canStart = canStartProcess(processToRun, result.finalStocks);
                
                // Check if any reusable resources are already in use
                if (canStart && reusableResources.find(processName) != reusableResources.end()) {
                    for (const auto& resource : reusableResources[processName]) {
                        auto inputIt = processToRun->inputs.find(resource);
                        if (inputIt != processToRun->inputs.end()) {
                            auto inUseIt = resourcesInUse.find(resource);
                            int inUseQty = (inUseIt != resourcesInUse.end()) ? inUseIt->second : 0;
                            int availableQty = result.finalStocks[resource] - inUseQty;
                            if (availableQty < inputIt->second) {
                                canStart = false;
                                break;
                            }
                        }
                    }
                }
                
                if (canStart && (currentTime + processToRun->nbCycle <= timeLimit)) {
                    // Consume inputs
                    for (const auto &[resource, quantity] : processToRun->inputs) {
                        result.finalStocks[resource] -= quantity;
                        
                        // Mark reusable resources as in use
                        if (reusableResources.find(processName) != reusableResources.end() &&
                            reusableResources[processName].find(resource) != reusableResources[processName].end()) {
                            resourcesInUse[resource] += quantity;
                        }
                    }

                    // Add to running processes queue
                    runningProcesses.push(
                        {processToRun,
                         currentTime + processToRun->nbCycle});
                    result.executionLog.push_back(
                        {currentTime, processName});
                    processStarted = true;
                    tryMoreProcesses = true;  // Try to start more processes in this cycle

                    // Move to next process in sequence only if successful
                    if (attemptedFromSequence) {
                        sequenceIndex++;
                    }
                } else {
                    // We couldn't start the process, only move to next if
                    // from sequence
                    if (!attemptedFromSequence) {
                        tryMoreProcesses = false;  // Stop trying opportunistic processes
                    } else {
                        // In sequence mode, move to next process if this one can't start
                        sequenceIndex++;
                        tryMoreProcesses = (sequenceIndex < processSequence.size());
                    }
                }
            } else {
                // Invalid process name, skip it
                if (attemptedFromSequence) {
                    sequenceIndex++;
                    tryMoreProcesses = (sequenceIndex < processSequence.size());
                }

                // Switch to opportunistic mode if sequence is exhausted
                if (!opportunisticMode &&
                    sequenceIndex >= processSequence.size()) {
                    opportunisticMode = true;
                    tryMoreProcesses = true;  // Try opportunistic mode
                }
            }
        }

        // Advance time if no process was started and no stocks were updated
        if (!processStarted && !stocksUpdated) {
            if (runningProcesses.empty()) {
                break;  // No more processes can run
            } else {
                // Jump to next process completion time
                int nextCompletionTime = runningProcesses.top().completionTime;
                currentTime = min(nextCompletionTime, timeLimit);
            }
        } else {
            // Process started or stocks updated, move to next cycle
            currentTime++;
        }

        if (currentTime > timeLimit) {
            currentTime = timeLimit;
            result.timeoutReached = true;
            break;
        }
    }

    // Process any remaining completions within the time limit
    while (!runningProcesses.empty() &&
           runningProcesses.top().completionTime <= timeLimit) {
        RunningProcess finishedProcess = runningProcesses.top();
        runningProcesses.pop();

        if (finishedProcess.processPtr) {
            for (const auto &[resource, quantity] :
                 finishedProcess.processPtr->outputs) {
                result.finalStocks[resource] += quantity;
            }
        }
    }

    // Calculate final cycle time
    if (!result.executionLog.empty()) {
        int maxEndTime = 0;
        for (const auto &[startTime, processName] : result.executionLog) {
            const Process *proc = getProcessByName(processName);
            if (proc) {
                maxEndTime = max(maxEndTime, startTime + proc->nbCycle);
            }
        }
        result.finalCycle = min(maxEndTime, timeLimit);
    }

    // Calculate fitness based on optimization goals
    result.fitness = 0.0;
    const auto &goals = config.getOptimizeGoal();

    // If no processes were executed, assign minimum fitness
    if (processesExecuted == 0 && result.executionLog.empty()) {
        result.fitness = -1.0e18;
        return result;
    }

    bool optimizeTime = false;
    for (const string &goal : goals) {
        if (goal == "time") {
            optimizeTime = true;
            result.fitness += 10000.0 / (1.0 + result.finalCycle);
        } else {
            // Add weighted value of specified resource to fitness
            double stockValue = 0;
            if (result.finalStocks.count(goal)) {
                stockValue = static_cast<double>(result.finalStocks.at(goal));
            }
            result.fitness += stockValue * 1000.0;
        }
    }

    // Small bonus for more processes executed (tie-breaker)
    result.fitness += processesExecuted;

    // Penalty if we're not optimizing time and we ran out of processes to run
    if (!optimizeTime && currentTime <= timeLimit && runningProcesses.empty() &&
        opportunisticMode) {
        // Less severe penalty
        if (result.fitness > 0) {
            result.fitness *= 0.9;
        }
    }

    return result;
}