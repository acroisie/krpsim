#include "ProcessManager.hpp"
#include "Config.hpp"
#include "Individual.hpp"
#include "Process.hpp"
#include <algorithm>
#include <iostream>
#include <limits>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <cmath>

using namespace std;

const int ProcessManager::POPULATION_SIZE = 200;
static const int MAX_SEQUENCE_LEN = 100;
static const int NUM_GENERATIONS = 100;
static const double MUTATION_RATE = 0.2;
static const int TOURNAMENT_SIZE = 3;

ProcessManager::ProcessManager(const Config &config, int delayLimit)
    : config_(config), delayLimit_(delayLimit), currentCycle_(0) {
    for (const auto &stock : config_.getStocks()) {
        currentStocks_[stock.name] = stock.quantity;
    }
}

// Used for testing if a sequence is realizable
bool ProcessManager::simulateSequence(const std::vector<std::string>& sequence, 
                                     std::map<std::string, int>& finalStocks,
                                     int& executedCount) {
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
        const string& procName = sequence[idx];
        const Process* proc = nullptr;
        
        // Find the process
        for (const auto& p : config_.getProcesses()) {
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
        for (const auto& input : proc->inputs) {
            if (stocks[input.first] < input.second) {
                canExecute = false;
                break;
            }
        }
        
        if (canExecute) {
            // Execute the process
            for (const auto& input : proc->inputs) {
                stocks[input.first] -= input.second;
            }
            runningProcesses.push_back({proc, cycle + proc->nbCycle});
            executedCount++;
            idx++;
        } else {
            // Can't execute, advance time if processes are running
            if (!runningProcesses.empty()) {
                int nextCompletion = delayLimit_;
                for (const auto& rp : runningProcesses) {
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
        for (const auto& rp : runningProcesses) {
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

void ProcessManager::createResourceDependencySequence(vector<string>& sequence, 
                                                     const vector<string>& allProcessNames) {
    random_device rd;
    mt19937 gen(rd());
    
    // Identify optimization goals
    const auto& goals = config_.getOptimizeGoal();
    
    // Build a proper dependency graph of processes
    // For each goal, backtrack through dependencies to find prerequisite processes
    
    // Map each resource to processes that produce it
    unordered_map<string, vector<const Process*>> resourceProducers;
    // Map each resource to processes that consume it
    unordered_map<string, vector<const Process*>> resourceConsumers;
    
    // Build the maps
    for (const auto& proc : config_.getProcesses()) {
        // Map outputs to this producer
        for (const auto& output : proc.outputs) {
            resourceProducers[output.first].push_back(&proc);
        }
        
        // Map inputs to this consumer
        for (const auto& input : proc.inputs) {
            resourceConsumers[input.first].push_back(&proc);
        }
    }
    
    // Identify terminal processes (those that produce goal resources)
    vector<const Process*> terminalProcesses;
    
    if (!goals.empty()) {
        for (const auto& goal : goals) {
            if (goal != "time") {
                // Find processes that directly produce this goal
                for (const auto& proc : config_.getProcesses()) {
                    for (const auto& output : proc.outputs) {
                        if (output.first == goal) {
                            terminalProcesses.push_back(&proc);
                            break;
                        }
                    }
                }
            }
        }
    }
    
    // If we have no terminal processes, consider all processes as potential terminals
    if (terminalProcesses.empty()) {
        for (const auto& proc : config_.getProcesses()) {
            terminalProcesses.push_back(&proc);
        }
    }
    
    // For each terminal process, find a chain of prerequisites
    vector<vector<const Process*>> processChains;
    
    for (const auto* terminal : terminalProcesses) {
        vector<const Process*> chain;
        chain.push_back(terminal);
        
        // Build a queue of resources needed by this chain
        queue<string> neededResources;
        unordered_set<string> processedResources;
        
        // Add all inputs of the terminal process
        for (const auto& input : terminal->inputs) {
            neededResources.push(input.first);
        }
        
        // Process the queue
        while (!neededResources.empty()) {
            string resource = neededResources.front();
            neededResources.pop();
            
            if (processedResources.count(resource)) continue;
            processedResources.insert(resource);
            
            // Find processes that produce this resource
            if (resourceProducers.count(resource)) {
                const auto& producers = resourceProducers[resource];
                
                if (!producers.empty()) {
                    // Pick a random producer for this resource
                    uniform_int_distribution<> prodDist(0, producers.size() - 1);
                    const Process* producer = producers[prodDist(gen)];
                    
                    // Add it to the chain
                    chain.push_back(producer);
                    
                    // Add its inputs to the needed resources
                    for (const auto& input : producer->inputs) {
                        if (!processedResources.count(input.first)) {
                            neededResources.push(input.first);
                        }
                    }
                }
            }
        }
        
        // Reverse the chain to get prerequisite -> goal order
        reverse(chain.begin(), chain.end());
        
        // Add to our set of chains
        processChains.push_back(chain);
    }
    
    // Now create a sequence by interleaving process chains
    // This ensures we build up to goals in a logical order
    
    // Determine how many times to repeat the chains
    uniform_int_distribution<> repeatDist(3, 6);
    int repetitions = repeatDist(gen);
    
    for (int i = 0; i < repetitions; i++) {
        // For each chain
        for (const auto& chain : processChains) {
            // Add each process in the chain
            for (const auto* proc : chain) {
                sequence.push_back(proc->name);
            }
        }
    }
    
    // Add some randomization
    uniform_int_distribution<> extraDist(0, 5);
    int extraProcesses = extraDist(gen);
    
    for (int i = 0; i < extraProcesses; i++) {
        uniform_int_distribution<> procDist(0, allProcessNames.size() - 1);
        sequence.push_back(allProcessNames[procDist(gen)]);
    }
}

void ProcessManager::initializePopulation() {
    population_.resize(POPULATION_SIZE);
    
    vector<string> allProcessNames;
    for (const auto &proc : config_.getProcesses()) {
        allProcessNames.push_back(proc.name);
    }

    random_device rd;
    mt19937 gen(rd());
    
    // Create initial population
    for (int i = 0; i < POPULATION_SIZE; i++) {
        vector<string> sequence;
        bool valid = false;
        int attempts = 0;
        
        // Make 5 attempts to create a valid sequence
        while (!valid && attempts < 5) {
            attempts++;
            sequence.clear();
            
            // Semi-intelligently build a sequence
            if (i < POPULATION_SIZE / 2) {
                // Create an intelligent sequence that follows resource dependencies
                createResourceDependencySequence(sequence, allProcessNames);
            } else {
                // Create a simple random sequence
                uniform_int_distribution<> lenDist(10, MAX_SEQUENCE_LEN);
                int seqLen = lenDist(gen);
                
                for (int j = 0; j < seqLen; j++) {
                    uniform_int_distribution<> procDist(0, (int)allProcessNames.size() - 1);
                    sequence.push_back(allProcessNames[procDist(gen)]);
                }
            }
            
            // Test if the sequence is at least partially executable
            map<string, int> finalStocks;
            int executedCount;
            valid = simulateSequence(sequence, finalStocks, executedCount) && executedCount > 0;
        }
        
        population_[i] = Individual(sequence);
    }
}

double ProcessManager::calculateFitness(Individual &individual) {
    // Reset state for simulation
    map<string, int> stocks;
    for (const auto &stock : config_.getStocks()) {
        stocks[stock.name] = stock.quantity;
    }
    executionLogs_.clear();
    currentCycle_ = 0;

    // Track initial stocks for profit calculation
    map<string, int> initialStocks = stocks;
    
    // Track resource usage to avoid excessive stockpiling
    map<string, int> resourcesUsed;
    map<string, int> resourcesProduced;
    
    struct RunningProcess {
        const Process *process;
        int completionCycle;
    };
    vector<RunningProcess> runningProcesses;

    // Count executed processes
    int executed = 0;
    size_t index = 0;
    
    while (currentCycle_ < delayLimit_ && index < individual.processSequence.size()) {
        // Complete finished processes
        auto it = runningProcesses.begin();
        while (it != runningProcesses.end()) {
            if (it->completionCycle <= currentCycle_) {
                for (const auto &out : it->process->outputs) {
                    stocks[out.first] += out.second;
                    resourcesProduced[out.first] += out.second;
                }
                it = runningProcesses.erase(it);
            } else {
                ++it;
            }
        }

        // Try to execute the current process
        const string &procName = individual.processSequence[index];
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
            if (stocks[input.first] < input.second) {
                canExecute = false;
                break;
            }
        }

        if (canExecute) {
            // Consume inputs
            for (const auto &input : chosen->inputs) {
                stocks[input.first] -= input.second;
                resourcesUsed[input.first] += input.second;
            }
            runningProcesses.push_back({chosen, currentCycle_ + chosen->nbCycle});
            executionLogs_.push_back({currentCycle_, chosen->name});
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
                currentCycle_ = nextCompletion;
            } else {
                // No running processes, skip this one
                index++;
            }
        }
    }

    // Complete any remaining processes
    while (!runningProcesses.empty() && currentCycle_ < delayLimit_) {
        int nextCompletion = delayLimit_;
        for (auto &rp : runningProcesses) {
            nextCompletion = min(nextCompletion, rp.completionCycle);
        }
        currentCycle_ = nextCompletion;

        auto it = runningProcesses.begin();
        while (it != runningProcesses.end()) {
            if (it->completionCycle <= currentCycle_) {
                for (const auto &out : it->process->outputs) {
                    stocks[out.first] += out.second;
                    resourcesProduced[out.first] += out.second;
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
        // For each goal
        for (const auto &goal : goals) {
            if (goal == "time") {
                // Minimize time - lower is better
                double timeFactor = 1.0 - ((double)currentCycle_ / delayLimit_);
                fitness += timeFactor * 100.0;
            } else {
                // Maximize goal stock quantity with proper aggressive normalization
                if (stocks.count(goal) > 0) {
                    // Very aggressive log normalization for values > 100
                    double normalizedValue = stocks[goal];
                    if (normalizedValue > 100) {
                        // Use log base 10 for really large values
                        normalizedValue = 100 + 20 * log10(normalizedValue / 100);
                    }
                    fitness += normalizedValue * 10.0;
                }
            }
        }
    }

    // Reward for resource utilization efficiency
    double utilizationScore = 0.0;
    for (const auto& used : resourcesUsed) {
        // Resources that were used in processes
        const string& resourceName = used.first;
        int amountUsed = used.second;
        // int amountProduced = resourcesProduced[resourceName];
        int finalAmount = stocks[resourceName];
        int initialAmount = initialStocks[resourceName];
        
        // Calculate how effectively this resource was used
        if (initialAmount > 0) {
            // For resources we started with - reward high usage percentage
            double usageRatio = (double)amountUsed / initialAmount;
            utilizationScore += min(1.0, usageRatio) * 50.0; // Cap at 100% usage
        }
        
        // Penalize excessive stockpiling of a resource
        if (finalAmount > initialAmount && resourceName != "dollars" && 
            resourceName != "reputation" && resourceName != "friendship") {
            // We're stockpiling this resource, which may be inefficient
            double stockpilingPenalty = 0.0;
            
            // Only penalize if it's not a goal resource
            bool isGoal = false;
            for (const auto& goal : goals) {
                if (goal == resourceName) {
                    isGoal = true;
                    break;
                }
            }
            
            if (!isGoal) {
                // Calculate penalty based on how much more than initial is stockpiled
                double excessRatio = (double)(finalAmount - initialAmount) / (initialAmount + 1);
                stockpilingPenalty = min(50.0, excessRatio * 10.0);
                utilizationScore -= stockpilingPenalty;
            }
        }
    }
    
    // Add utilization score to fitness
    fitness += utilizationScore;
    
    // Analyze process chain completions
    // Determine if key process chains were completed (e.g., buy_boat -> find_treasure -> sell_treasure)
    unordered_map<string, bool> processExecuted;
    for (const auto& log : executionLogs_) {
        processExecuted[log.second] = true;
    }
    
    // Check for specific process chains in a generic way
    // Look for complete chains where each process creates inputs for the next
    double chainCompletionBonus = 0.0;
    
    // For each process that produces a goal resource
    for (const auto &goal : goals) {
        if (goal != "time") {
            for (const auto &proc : config_.getProcesses()) {
                bool producesGoal = false;
                for (const auto &out : proc.outputs) {
                    if (out.first == goal) {
                        producesGoal = true;
                        break;
                    }
                }
                
                if (producesGoal && processExecuted[proc.name]) {
                    // This process produces the goal and was executed
                    // Find its prerequisites and check if they were executed
                    vector<string> prerequisites;
                    for (const auto &input : proc.inputs) {
                        // Find processes that produce this input
                        for (const auto &otherProc : config_.getProcesses()) {
                            for (const auto &out : otherProc.outputs) {
                                if (out.first == input.first) {
                                    prerequisites.push_back(otherProc.name);
                                }
                            }
                        }
                    }
                    
                    // Check if all prerequisites were executed
                    bool allPrereqsExecuted = true;
                    for (const auto &prereq : prerequisites) {
                        if (!processExecuted[prereq]) {
                            allPrereqsExecuted = false;
                            break;
                        }
                    }
                    
                    // Bonus for complete process chains
                    if (allPrereqsExecuted && !prerequisites.empty()) {
                        chainCompletionBonus += 200.0;
                    }
                }
            }
        }
    }
    
    // Add chain completion bonus
    fitness += chainCompletionBonus;
    
    // Reward for number of processes executed
    fitness += executed * 10.0;
    
    return fitness;
}

vector<const Process*> ProcessManager::getRunnableProcesses() {
    vector<const Process*> result;
    for (const auto &proc : config_.getProcesses()) {
        bool canRun = true;
        for (const auto &in : proc.inputs) {
            if (currentStocks_.count(in.first) == 0 || currentStocks_[in.first] < in.second) {
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

bool ProcessManager::executeProcess(const Process *process) {
    if (!process) return false;
    
    // Verify we have enough resources
    for (const auto &in : process->inputs) {
        if (currentStocks_.count(in.first) == 0 || currentStocks_[in.first] < in.second) {
            return false;
        }
    }
    
    // Consume inputs
    for (const auto &in : process->inputs) {
        currentStocks_[in.first] -= in.second;
    }
    
    // Produce outputs
    for (const auto &out : process->outputs) {
        currentStocks_[out.first] += out.second;
    }
    
    return true;
}

void ProcessManager::evaluateFitness() {
    for (auto &indiv : population_) {
        indiv.fitness = calculateFitness(indiv);
    }
}

vector<Individual> ProcessManager::selectParents(const vector<Individual> &population) {
    vector<Individual> parents;
    if (population.empty()) {
        return parents;
    }

    random_device rd;
    mt19937 gen(rd());
    
    int parentsToSelect = min(POPULATION_SIZE, (int)population.size());
    
    // Elitism: always keep the best individual
    auto bestIndiv = findBestIndividual(population);
    parents.push_back(bestIndiv);
    
    // Select the rest using tournament selection
    for (int i = 1; i < parentsToSelect; ++i) {
        // Select TOURNAMENT_SIZE random individuals
        vector<const Individual*> tournament;
        for (int j = 0; j < TOURNAMENT_SIZE; ++j) {
            uniform_int_distribution<> dist(0, (int)population.size() - 1);
            tournament.push_back(&population[dist(gen)]);
        }
        
        // Find the best in tournament
        const Individual* best = tournament[0];
        for (int j = 1; j < TOURNAMENT_SIZE; ++j) {
            if (tournament[j]->fitness > best->fitness) {
                best = tournament[j];
            }
        }
        
        parents.push_back(*best);
    }
    
    return parents;
}

pair<Individual, Individual> ProcessManager::crossover(const Individual &parent1, const Individual &parent2) {
    if (parent1.processSequence.empty() || parent2.processSequence.empty()) {
        return make_pair(parent1, parent2);
    }
    
    random_device rd;
    mt19937 gen(rd());
    
    size_t minSize = min(parent1.processSequence.size(), parent2.processSequence.size());
    if (minSize < 2) {
        return make_pair(parent1, parent2);
    }

    uniform_int_distribution<> pointDist(1, (int)minSize - 1);
    size_t crossoverPoint = pointDist(gen);
    
    vector<string> child1, child2;
    
    // First child gets first part of parent1, second part of parent2
    child1.insert(child1.end(), 
                 parent1.processSequence.begin(), 
                 parent1.processSequence.begin() + crossoverPoint);
    
    if (crossoverPoint < parent2.processSequence.size()) {
        child1.insert(child1.end(), 
                    parent2.processSequence.begin() + crossoverPoint, 
                    parent2.processSequence.end());
    }
    
    // Second child gets first part of parent2, second part of parent1
    child2.insert(child2.end(), 
                 parent2.processSequence.begin(), 
                 parent2.processSequence.begin() + crossoverPoint);
    
    if (crossoverPoint < parent1.processSequence.size()) {
        child2.insert(child2.end(), 
                    parent1.processSequence.begin() + crossoverPoint, 
                    parent1.processSequence.end());
    }

    return make_pair(Individual(child1), Individual(child2));
}

Individual ProcessManager::mutate(const Individual &individual, double mutationRate) {
    auto allProcs = config_.getProcesses();
    if (allProcs.empty() || individual.processSequence.empty()) {
        return individual;
    }

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> distProb(0.0, 1.0);
    
    vector<string> mutatedSeq = individual.processSequence;
    
    // Collect process names for convenience
    vector<string> allProcessNames;
    for (const auto &proc : allProcs) {
        allProcessNames.push_back(proc.name);
    }
    
    // Build maps of resource dependencies for smarter mutations
    unordered_map<string, vector<string>> resourceProducers;
    unordered_map<string, vector<string>> resourceConsumers;
    
    for (const auto& proc : allProcs) {
        // Record outputs (what this process produces)
        for (const auto& output : proc.outputs) {
            resourceProducers[output.first].push_back(proc.name);
        }
        
        // Record inputs (what this process consumes)
        for (const auto& input : proc.inputs) {
            resourceConsumers[input.first].push_back(proc.name);
        }
    }
    
    // Point mutations - different strategies based on position in sequence
    for (size_t i = 0; i < mutatedSeq.size(); ++i) {
        if (distProb(gen) < mutationRate) {
            // 50% chance for regular random mutation
            if (distProb(gen) < 0.5) {
                uniform_int_distribution<> procDist(0, (int)allProcessNames.size() - 1);
                mutatedSeq[i] = allProcessNames[procDist(gen)];
            }
            // 50% chance for intelligent mutation
            else {
                const string& currentProc = mutatedSeq[i];
                
                // Find the current process
                const Process* current = nullptr;
                for (const auto& proc : allProcs) {
                    if (proc.name == currentProc) {
                        current = &proc;
                        break;
                    }
                }
                
                if (current) {
                    // Choose a mutation strategy based on process position
                    
                    // 1. Find processes that produce same outputs
                    vector<string> similarProducers;
                    for (const auto& output : current->outputs) {
                        for (const string& procName : resourceProducers[output.first]) {
                            if (procName != currentProc) {
                                similarProducers.push_back(procName);
                            }
                        }
                    }
                    
                    // 2. Find processes that consume same inputs
                    vector<string> similarConsumers;
                    for (const auto& input : current->inputs) {
                        for (const string& procName : resourceConsumers[input.first]) {
                            if (procName != currentProc) {
                                similarConsumers.push_back(procName);
                            }
                        }
                    }
                    
                    // 3. Find processes that produce this process's inputs
                    vector<string> inputProducers;
                    for (const auto& input : current->inputs) {
                        for (const string& procName : resourceProducers[input.first]) {
                            inputProducers.push_back(procName);
                        }
                    }
                    
                    // 4. Find processes that consume this process's outputs
                    vector<string> outputConsumers;
                    for (const auto& output : current->outputs) {
                        for (const string& procName : resourceConsumers[output.first]) {
                            outputConsumers.push_back(procName);
                        }
                    }
                    
                    // Choose one of these lists based on process position in sequence
                    vector<string>* candidateList = nullptr;
                    
                    if (i < mutatedSeq.size() / 3) {
                        // Early in sequence - prefer input producers
                        if (!inputProducers.empty()) {
                            candidateList = &inputProducers;
                        } else if (!similarProducers.empty()) {
                            candidateList = &similarProducers;
                        }
                    } else if (i > 2 * mutatedSeq.size() / 3) {
                        // Late in sequence - prefer output consumers
                        if (!outputConsumers.empty()) {
                            candidateList = &outputConsumers;
                        } else if (!similarConsumers.empty()) {
                            candidateList = &similarConsumers;
                        }
                    } else {
                        // Middle of sequence - mix it up
                        uniform_int_distribution<> typeDist(0, 3);
                        int type = typeDist(gen);
                        if (type == 0 && !inputProducers.empty()) {
                            candidateList = &inputProducers;
                        } else if (type == 1 && !outputConsumers.empty()) {
                            candidateList = &outputConsumers;
                        } else if (type == 2 && !similarProducers.empty()) {
                            candidateList = &similarProducers;
                        } else if (type == 3 && !similarConsumers.empty()) {
                            candidateList = &similarConsumers;
                        }
                    }
                    
                    // If we found a suitable list, pick a process from it
                    if (candidateList && !candidateList->empty()) {
                        uniform_int_distribution<> candidateDist(0, candidateList->size() - 1);
                        mutatedSeq[i] = (*candidateList)[candidateDist(gen)];
                    } else {
                        // Fall back to random mutation
                        uniform_int_distribution<> procDist(0, (int)allProcessNames.size() - 1);
                        mutatedSeq[i] = allProcessNames[procDist(gen)];
                    }
                } else {
                    // Current process not found, use random mutation
                    uniform_int_distribution<> procDist(0, (int)allProcessNames.size() - 1);
                    mutatedSeq[i] = allProcessNames[procDist(gen)];
                }
            }
        }
    }
    
    // Structural mutation with more variety
    if (distProb(gen) < 0.3) {
        uniform_int_distribution<> mutType(0, 4);  // Added new mutation types
        int mutation = mutType(gen);
        
        switch (mutation) {
            case 0: // Insert
                if (mutatedSeq.size() < MAX_SEQUENCE_LEN) {
                    uniform_int_distribution<> posDist(0, (int)mutatedSeq.size());
                    uniform_int_distribution<> procDist(0, (int)allProcessNames.size() - 1);
                    mutatedSeq.insert(mutatedSeq.begin() + posDist(gen), 
                                     allProcessNames[procDist(gen)]);
                }
                break;
                
            case 1: // Delete
                if (mutatedSeq.size() > 1) {
                    uniform_int_distribution<> posDist(0, (int)mutatedSeq.size() - 1);
                    mutatedSeq.erase(mutatedSeq.begin() + posDist(gen));
                }
                break;
                
            case 2: // Swap
                if (mutatedSeq.size() >= 2) {
                    uniform_int_distribution<> posDist(0, (int)mutatedSeq.size() - 1);
                    size_t pos1 = posDist(gen);
                    size_t pos2 = posDist(gen);
                    swap(mutatedSeq[pos1], mutatedSeq[pos2]);
                }
                break;
                
            case 3: // Duplicate segment
                if (mutatedSeq.size() >= 3 && mutatedSeq.size() < MAX_SEQUENCE_LEN - 3) {
                    // Pick a small segment to duplicate
                    uniform_int_distribution<> startDist(0, (int)mutatedSeq.size() - 3);
                    uniform_int_distribution<> lenDist(1, 3);
                    
                    int start = startDist(gen);
                    int len = min(lenDist(gen), (int)mutatedSeq.size() - start);
                    
                    // Create the segment to duplicate
                    vector<string> segment(mutatedSeq.begin() + start, 
                                          mutatedSeq.begin() + start + len);
                    
                    // Insert it at a random position
                    uniform_int_distribution<> insDist(0, (int)mutatedSeq.size());
                    int insertPos = insDist(gen);
                    
                    mutatedSeq.insert(mutatedSeq.begin() + insertPos, 
                                     segment.begin(), segment.end());
                }
                break;
                
            case 4: // Insert process chain
                if (mutatedSeq.size() < MAX_SEQUENCE_LEN - 5) {
                    // Create a small chain of related processes
                    vector<string> chain;
                    createResourceDependencySequence(chain, allProcessNames);
                    
                    // Limit the chain size
                    if (chain.size() > 5) {
                        chain.resize(5);
                    }
                    
                    // Insert the chain at a random position
                    if (!chain.empty()) {
                        uniform_int_distribution<> insDist(0, (int)mutatedSeq.size());
                        int insertPos = insDist(gen);
                        
                        mutatedSeq.insert(mutatedSeq.begin() + insertPos, 
                                         chain.begin(), chain.end());
                    }
                }
                break;
        }
    }
    
    return Individual(mutatedSeq);
}

void ProcessManager::createDiversityRecoveryPopulation(vector<Individual>& newPopulation, 
                                                     const vector<string>& allProcessNames) {
    random_device rd;
    mt19937 gen(rd());
    
    // Add individuals with more sophisticated diversity strategies
    
    // Strategy 1: Process chain focused individuals
    int chainFocusedCount = POPULATION_SIZE / 6;
    for (int j = 0; j < chainFocusedCount && newPopulation.size() < POPULATION_SIZE; j++) {
        vector<string> chainSeq;
        createResourceDependencySequence(chainSeq, allProcessNames);
        
        // Test if the sequence is at least partially executable
        map<string, int> finalStocks;
        int executedCount;
        if (simulateSequence(chainSeq, finalStocks, executedCount) && executedCount > 0) {
            newPopulation.push_back(Individual(chainSeq));
        }
    }
    
    // Strategy 2: Greedy resource acquisition individuals
    // (Focused on acquiring as many resources as possible first)
    int greedyCount = POPULATION_SIZE / 6;
    for (int j = 0; j < greedyCount && newPopulation.size() < POPULATION_SIZE; j++) {
        vector<string> greedySeq;
        
        // Find processes that are executable from the initial state
        vector<const Process*> initialProcesses;
        map<string, int> initialStocks;
        for (const auto& stock : config_.getStocks()) {
            initialStocks[stock.name] = stock.quantity;
        }
        
        for (const auto& proc : config_.getProcesses()) {
            bool canExecute = true;
            for (const auto& input : proc.inputs) {
                if (initialStocks[input.first] < input.second) {
                    canExecute = false;
                    break;
                }
            }
            if (canExecute) {
                initialProcesses.push_back(&proc);
            }
        }
        
        // Add these processes multiple times to focus on resource acquisition
        if (!initialProcesses.empty()) {
            for (int k = 0; k < 20; k++) {
                uniform_int_distribution<> procDist(0, initialProcesses.size() - 1);
                greedySeq.push_back(initialProcesses[procDist(gen)]->name);
            }
            
            // Then add some random processes to use the resources
            for (int k = 0; k < 20; k++) {
                uniform_int_distribution<> procDist(0, allProcessNames.size() - 1);
                greedySeq.push_back(allProcessNames[procDist(gen)]);
            }
            
            // Test if the sequence is executable
            map<string, int> finalStocks;
            int executedCount;
            if (simulateSequence(greedySeq, finalStocks, executedCount) && executedCount > 0) {
                newPopulation.push_back(Individual(greedySeq));
            }
        }
    }
    
    // Strategy 3: Goal-focused individuals
    // (Directly targetting goal resources)
    int goalFocusedCount = POPULATION_SIZE / 6;
    const auto& goals = config_.getOptimizeGoal();
    if (!goals.empty()) {
        for (int j = 0; j < goalFocusedCount && newPopulation.size() < POPULATION_SIZE; j++) {
            vector<string> goalSeq;
            
            // Find processes that directly produce goal resources
            vector<const Process*> goalProducers;
            for (const auto& goal : goals) {
                if (goal != "time") {
                    for (const auto& proc : config_.getProcesses()) {
                        for (const auto& output : proc.outputs) {
                            if (output.first == goal) {
                                goalProducers.push_back(&proc);
                                break;
                            }
                        }
                    }
                }
            }
            
            // If we have goal producers, create sequences focused on them
            if (!goalProducers.empty()) {
                // First add processes that might produce prerequisites
                for (int k = 0; k < 15; k++) {
                    uniform_int_distribution<> procDist(0, allProcessNames.size() - 1);
                    goalSeq.push_back(allProcessNames[procDist(gen)]);
                }
                
                // Then heavily focus on the goal-producing processes
                for (int k = 0; k < 30; k++) {
                    uniform_int_distribution<> procDist(0, goalProducers.size() - 1);
                    goalSeq.push_back(goalProducers[procDist(gen)]->name);
                }
                
                // Test if the sequence is executable
                map<string, int> finalStocks;
                int executedCount;
                if (simulateSequence(goalSeq, finalStocks, executedCount) && executedCount > 0) {
                    newPopulation.push_back(Individual(goalSeq));
                }
            }
        }
    }
    
    // Strategy 4: Completely random individuals with validation
    int randomCount = POPULATION_SIZE / 6;
    for (int j = 0; j < randomCount && newPopulation.size() < POPULATION_SIZE; j++) {
        vector<string> randomSeq;
        uniform_int_distribution<> lenDist(10, MAX_SEQUENCE_LEN);
        int seqLen = lenDist(gen);
        
        for (int k = 0; k < seqLen; k++) {
            uniform_int_distribution<> procDist(0, allProcessNames.size() - 1);
            randomSeq.push_back(allProcessNames[procDist(gen)]);
        }
        
        // Test if the sequence is executable
        map<string, int> finalStocks;
        int executedCount;
        if (simulateSequence(randomSeq, finalStocks, executedCount) && executedCount > 0) {
            newPopulation.push_back(Individual(randomSeq));
        }
    }
    
    // Strategy 5: Repeat & rotate proven successful patterns
    // (Take best solutions and create variations)
    if (!population_.empty() && newPopulation.size() < POPULATION_SIZE) {
        vector<Individual> sortedPopulation = population_;
        sort(sortedPopulation.begin(), sortedPopulation.end(), 
             [](const Individual& a, const Individual& b) { 
                 return a.fitness > b.fitness; 
             });
             
        int bestCount = min(5, (int)sortedPopulation.size());
        int variationsPerBest = (POPULATION_SIZE - newPopulation.size()) / (bestCount + 1);
        
        for (int j = 0; j < bestCount; j++) {
            const auto& bestIndiv = sortedPopulation[j];
            
            for (int k = 0; k < variationsPerBest && newPopulation.size() < POPULATION_SIZE; k++) {
                // Create variations by shuffling segments of the sequence
                vector<string> variantSeq = bestIndiv.processSequence;
                
                if (variantSeq.size() >= 10) {
                    // Choose a random segment to shuffle
                    uniform_int_distribution<> startDist(0, variantSeq.size() - 10);
                    uniform_int_distribution<> lengthDist(5, 10);
                    
                    int start = startDist(gen);
                    int length = min(lengthDist(gen), (int)variantSeq.size() - start);
                    
                    // Shuffle the segment
                    shuffle(variantSeq.begin() + start, variantSeq.begin() + start + length, gen);
                    
                    // Test if the modified sequence is still executable
                    map<string, int> finalStocks;
                    int executedCount;
                    if (simulateSequence(variantSeq, finalStocks, executedCount) && executedCount > 0) {
                        newPopulation.push_back(Individual(variantSeq));
                    }
                }
            }
        }
    }
}

Individual ProcessManager::findBestIndividual(const vector<Individual> &population) {
    if (population.empty()) {
        return Individual();
    }
    
    Individual best = population.front();
    for (const auto &ind : population) {
        if (ind.fitness > best.fitness) {
            best = ind;
        }
    }
    return best;
}

void ProcessManager::runGeneticAlgorithm() {
    cout << "Starting genetic algorithm optimization..." << endl;
    initializePopulation();
    
    evaluateFitness();
    
    // Track the best solution found
    if (!population_.empty()) {
        bestSolution_ = findBestIndividual(population_);
    }

    // Track previous best fitness to detect stagnation
    double prevBestFitness = bestSolution_.fitness;
    int stagnationCounter = 0;
    
    // Collect all process names for convenience
    vector<string> allProcessNames;
    for (const auto &proc : config_.getProcesses()) {
        allProcessNames.push_back(proc.name);
    }
    
    // Run for specified number of generations
    for (int generation = 0; generation < NUM_GENERATIONS; ++generation) {
        // Select parents for next generation
        vector<Individual> parents = selectParents(population_);
        vector<Individual> newPopulation;

        // Elitism - keep the best solution
        if (!population_.empty()) {
            auto best = findBestIndividual(population_);
            if (best.fitness > bestSolution_.fitness) {
                bestSolution_ = best;
            }
            newPopulation.push_back(best);
            
            // Also keep a few more top solutions
            vector<Individual> sortedPopulation = population_;
            sort(sortedPopulation.begin(), sortedPopulation.end(), 
                 [](const Individual& a, const Individual& b) { 
                     return a.fitness > b.fitness; 
                 });
                 
            int eliteCount = min(5, (int)sortedPopulation.size());
            for (int j = 1; j < eliteCount; j++) {
                newPopulation.push_back(sortedPopulation[j]);
            }
        }
        
        // Stagnation detection and recovery
        double currentBestFitness = findBestIndividual(population_).fitness;
        if (currentBestFitness <= prevBestFitness) {
            stagnationCounter++;
            cout << "No improvement for " << stagnationCounter << " generations (fitness: " << currentBestFitness << ")" << endl;
        } else {
            cout << "Fitness improved from " << prevBestFitness << " to " << currentBestFitness << endl;
            stagnationCounter = 0;
            prevBestFitness = currentBestFitness;
        }
        
        // If stagnating for 5 generations, inject more diversity
        if (stagnationCounter >= 5) {
            cout << "Stagnation detected, injecting diversity in generation " << generation << endl;
            stagnationCounter = 0;
            
            // Add diversity recovery population
            createDiversityRecoveryPopulation(newPopulation, allProcessNames);
        }

        // Create new population through crossover and mutation
        while (newPopulation.size() < population_.size() && parents.size() >= 2) {
            // Select two parents randomly
            random_device rd;
            mt19937 gen(rd());
            uniform_int_distribution<> dist(0, (int)parents.size() - 1);
            
            size_t idx1 = dist(gen);
            size_t idx2 = dist(gen);
            while (idx2 == idx1 && parents.size() > 1) {
                idx2 = dist(gen);
            }

            auto children = crossover(parents[idx1], parents[idx2]);
            auto child1 = mutate(children.first, MUTATION_RATE);
            auto child2 = mutate(children.second, MUTATION_RATE);

            if (newPopulation.size() < population_.size()) {
                // Validate children before adding them
                map<string, int> finalStocks;
                int executedCount;
                if (simulateSequence(child1.processSequence, finalStocks, executedCount) && executedCount > 0) {
                    newPopulation.push_back(child1);
                }
            }
            
            if (newPopulation.size() < population_.size()) {
                map<string, int> finalStocks;
                int executedCount;
                if (simulateSequence(child2.processSequence, finalStocks, executedCount) && executedCount > 0) {
                    newPopulation.push_back(child2);
                }
            }
        }

        // If we don't have enough individuals, add random ones
        while (newPopulation.size() < population_.size()) {
            random_device rd;
            mt19937 gen(rd());
            
            vector<string> randomSeq;
            uniform_int_distribution<> lenDist(10, MAX_SEQUENCE_LEN);
            int seqLen = lenDist(gen);
            
            for (int j = 0; j < seqLen; j++) {
                uniform_int_distribution<> procDist(0, (int)allProcessNames.size() - 1);
                randomSeq.push_back(allProcessNames[procDist(gen)]);
            }
            
            // Test if the sequence is at least partially executable
            map<string, int> finalStocks;
            int executedCount;
            if (simulateSequence(randomSeq, finalStocks, executedCount) && executedCount > 0) {
                newPopulation.push_back(Individual(randomSeq));
            }
        }

        // Replace population
        if (!newPopulation.empty()) {
            population_ = newPopulation;
        }
        
        // Evaluate fitness of new population
        evaluateFitness();
        
        // Track the best solution in this generation
        auto currentBest = findBestIndividual(population_);
        
        // Update best solution if improvement found
        if (currentBest.fitness > bestSolution_.fitness) {
            bestSolution_ = currentBest;
        }
        
        // Output progress
        cout << "Generation " << generation + 1 << ": Best fitness = " << bestSolution_.fitness << endl;
    }

    cout << "Genetic algorithm completed. Generating final solution." << endl;
    generateOutput();
}

void ProcessManager::generateOutput() {
    // Reset state for final simulation
    currentStocks_.clear();
    for (const auto &stock : config_.getStocks()) {
        currentStocks_[stock.name] = stock.quantity;
    }
    executionLogs_.clear();
    currentCycle_ = 0;

    struct RunningProcess {
        const Process *process;
        int completionCycle;
    };
    vector<RunningProcess> runningProcesses;

    const auto &sequence = bestSolution_.processSequence;
    size_t idx = 0;

    while (currentCycle_ < delayLimit_ && idx < sequence.size()) {
        // Complete finished processes
        auto it = runningProcesses.begin();
        while (it != runningProcesses.end()) {
            if (it->completionCycle <= currentCycle_) {
                for (auto &out : it->process->outputs) {
                    currentStocks_[out.first] += out.second;
                }
                it = runningProcesses.erase(it);
            } else {
                ++it;
            }
        }

        // Try to run next process
        const string &procName = sequence[idx];
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
            idx++;
            continue;
        }

        // Check if we can execute it
        bool canExecute = true;
        for (const auto &input : chosen->inputs) {
            if (currentStocks_[input.first] < input.second) {
                canExecute = false;
                break;
            }
        }

        if (canExecute) {
            // Consume inputs
            for (auto &in : chosen->inputs) {
                currentStocks_[in.first] -= in.second;
            }
            runningProcesses.push_back({chosen, currentCycle_ + chosen->nbCycle});
            executionLogs_.push_back({currentCycle_, chosen->name});
            ++idx;
        } else {
            if (!runningProcesses.empty()) {
                // Fast-forward to next process completion
                int nextCompletion = delayLimit_;
                for (auto &rp : runningProcesses) {
                    nextCompletion = min(nextCompletion, rp.completionCycle);
                }
                currentCycle_ = nextCompletion;
            } else {
                // Skip this process and move to next
                ++idx;
                ++currentCycle_;
            }
        }
    }

    // Complete any remaining processes
    while (!runningProcesses.empty() && currentCycle_ < delayLimit_) {
        int nextCompletion = delayLimit_;
        for (auto &rp : runningProcesses) {
            nextCompletion = min(nextCompletion, rp.completionCycle);
        }
        currentCycle_ = nextCompletion;

        auto it = runningProcesses.begin();
        while (it != runningProcesses.end()) {
            if (it->completionCycle <= currentCycle_) {
                for (auto &out : it->process->outputs) {
                    currentStocks_[out.first] += out.second;
                }
                it = runningProcesses.erase(it);
            } else {
                ++it;
            }
        }
    }

    // Output results
    cout << "Nice file! " << config_.getProcesses().size() << " processes, " << config_.getStocks().size()
         << " stocks, " << config_.getOptimizeGoal().size() << " to optimize" << endl;

    cout << "Main walk:" << endl;
    for (auto &log : executionLogs_) {
        cout << log.first << ":" << log.second << endl;
    }

    cout << "No more process doable at time " << currentCycle_ << endl;

    cout << "Stock:" << endl;
    for (auto &stock : currentStocks_) {
        cout << stock.first << " => " << stock.second << endl;
    }
}