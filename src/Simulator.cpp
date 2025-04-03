#include "Simulator.hpp"
#include <algorithm>

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
                for (const auto &[resource, quantity] : finishedProcess.processPtr->inputs) {
                    auto outputIt = finishedProcess.processPtr->outputs.find(resource);
                    if (outputIt != finishedProcess.processPtr->outputs.end() && outputIt->second == quantity) {
                        // This is a reusable resource (like the pan)
                        resourcesInUse[resource] -= quantity;
                        if (resourcesInUse[resource] <= 0) {
                            resourcesInUse.erase(resource);
                        }
                    }
                }
                
                stocksUpdated = true;
            }
        }

        if (currentTime >= timeLimit) {
            result.timeoutReached = true;
            break;
        }

        bool processStarted = false;
        if (currentTime < timeLimit) {
            bool tryMore = true;
            while (tryMore) {
                tryMore = false;
                const Process *processToRun = nullptr;
                string processName;
                bool attemptedFromSequence = false;

                // Try to get a process from the sequence or opportunistically
                if (!opportunisticMode &&
                    sequenceIndex < processSequence.size()) {
                    processName = processSequence[sequenceIndex];
                    processToRun = getProcessByName(processName);
                    attemptedFromSequence = true;
                } else {
                    opportunisticMode = true;
                    for (const auto &[name, process] : processMap) {
                        if (canStartProcess(process, result.finalStocks)) {
                            processToRun = process;
                            processName = name;
                            break;
                        }
                    }
                }

                if (processToRun) {
                    // Check if resources are available and not in use
                    bool canStart = canStartProcess(processToRun, result.finalStocks);
                    
                    // Check if any reusable resources are already in use
                    if (canStart) {
                        for (const auto &[resource, quantity] : processToRun->inputs) {
                            auto outputIt = processToRun->outputs.find(resource);
                            if (outputIt != processToRun->outputs.end() && outputIt->second == quantity) {
                                // This is a reusable resource (like the pan)
                                auto inUseIt = resourcesInUse.find(resource);
                                if (inUseIt != resourcesInUse.end() && inUseIt->second > 0) {
                                    // This resource is already in use
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
                            auto outputIt = processToRun->outputs.find(resource);
                            if (outputIt != processToRun->outputs.end() && outputIt->second == quantity) {
                                resourcesInUse[resource] += quantity;
                            }
                        }

                        // Add to running processes queue
                        runningProcesses.push(
                            {processToRun,
                             currentTime + processToRun->nbCycle});
                        result.executionLog.push_back(
                            {currentTime, processName});
                        processesExecuted++;
                        processStarted = true;
                        tryMore = true;

                        // Move to next process in sequence only if successful
                        if (attemptedFromSequence) {
                            sequenceIndex++;
                        }
                    } else {
                        // We couldn't start the process, only move to next if
                        // opportunistic
                        if (!attemptedFromSequence) {
                            tryMore = false;
                        } else {
                            // In sequence mode, move to next process if this one can't start
                            sequenceIndex++;
                            tryMore = (sequenceIndex < processSequence.size());
                        }
                    }
                } else {
                    // Invalid process name, skip it
                    if (attemptedFromSequence) {
                        sequenceIndex++;
                        tryMore = (sequenceIndex < processSequence.size());
                    }

                    // Switch to opportunistic mode if sequence is exhausted
                    if (!opportunisticMode &&
                        sequenceIndex >= processSequence.size()) {
                        opportunisticMode = true;
                        tryMore = true;
                    }
                }
            }
        }

        // Advance time if no process was started and no stocks were updated
        if (!processStarted && !stocksUpdated) {
            if (runningProcesses.empty()) {
                break;
            } else {
                int nextCompletionTime = runningProcesses.top().completionTime;
                if (nextCompletionTime > timeLimit) {
                    currentTime = timeLimit;
                    result.timeoutReached = true;
                } else {
                    // Make sure time advances at least by 1 to avoid infinite
                    // loops
                    currentTime = max(currentTime + 1, nextCompletionTime);
                }
            }
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
    result.fitness += processesExecuted * 1e-9;

    // Penalty if we're not optimizing time and we ran out of processes to run
    if (!optimizeTime && currentTime <= timeLimit && runningProcesses.empty() &&
        opportunisticMode) {
        result.fitness *= 0.5;
    }

    return result;
}