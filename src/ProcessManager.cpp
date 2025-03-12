#include "ProcessManager.hpp"
#include "Config.hpp"
#include "Individual.hpp"
#include "Process.hpp"
#include <iostream>
#include <limits>
#include <vector>
#include <string>
#include <algorithm>
#include <random>

using namespace std;

const int ProcessManager::POPULATION_SIZE = 20;

ProcessManager::ProcessManager(const Config &config, int delayLimit)
    : config_(config), delayLimit_(delayLimit), currentCycle_(0) {
    for (const auto &stock : config_.getStocks()) {
        currentStocks_[stock.name] = stock.quantity;
    }
}

void ProcessManager::initializePopulation() {
    population_.resize(POPULATION_SIZE);
    vector<string> availableProcessNames;
    for (const auto &process : config_.getProcesses()) {
        availableProcessNames.push_back(process.name);
    }
    random_device rd;
    mt19937 generator(rd());
    uniform_int_distribution<> lengthDistribution(1, delayLimit_);
    uniform_int_distribution<> processDistribution(0, availableProcessNames.size() - 1);
    for (int i = 0; i < POPULATION_SIZE; ++i) {
        vector<string> randomSequence;
        int sequenceLength = lengthDistribution(generator);
        for (int j = 0; j < sequenceLength; ++j) {
            if (!availableProcessNames.empty()) {
                int processIndex = processDistribution(generator);
                randomSequence.push_back(availableProcessNames[processIndex]);
            }
        }
        population_[i] = Individual(randomSequence);
    }
}

void ProcessManager::evaluateFitness() {
    for (int i = 0; i < POPULATION_SIZE; ++i) {
        population_[i].fitness = calculateFitness(population_[i]);
    }
}

vector<Individual> ProcessManager::selection() {
    vector<Individual> parents = selectParents(population_);
    cout << "Parents selected." << endl;
    return parents;
}

Individual ProcessManager::findBestIndividual(const std::vector<Individual>& population) {
    if (population.empty()) {
        return Individual(std::vector<std::string>());
    }
    
    Individual best = population[0];
    for (const auto& individual : population) {
        if (individual.fitness > best.fitness) {
            best = individual;
        }
    }
    
    return best;
}

void ProcessManager::runGeneticAlgorithm() {
    std::cout << "Running genetic algorithm... (runGeneticAlgorithm() START)" << std::endl;
    
    // Initialize population with random process sequences
    initializePopulation();
    std::cout << "Initial population generated." << std::endl;
    
    // Set parameters for the genetic algorithm
    const int NUM_GENERATIONS = 20; // Increase from 5 to 20
    const double MUTATION_RATE = 0.1;
    
    std::cout << "Starting evolution over " << NUM_GENERATIONS << " generations..." << std::endl;
    
    // Evaluate fitness of the initial population
    evaluateFitness();
    
    for (int generation = 0; generation < NUM_GENERATIONS; ++generation) {
        std::cout << "--- Generation " << generation << " ---" << std::endl;
        
        // Select parents based on fitness
        std::vector<Individual> parents = selectParents(population_);
        std::cout << "Parents selected." << std::endl;
        
        // Create new population through crossover and mutation
        std::vector<Individual> newPopulation;
        
        // Elitism: Keep the best individual
        if (!population_.empty()) {
            Individual bestIndividual = findBestIndividual(population_);
            newPopulation.push_back(bestIndividual);
        }
        
        // Generate children until we fill the population
        while (newPopulation.size() < population_.size() && parents.size() >= 2) {
            // Select two parents
            size_t idx1 = rand() % parents.size();
            size_t idx2 = rand() % parents.size();
            while (idx2 == idx1 && parents.size() > 1) {
                idx2 = rand() % parents.size();
            }
            
            // Perform crossover
            auto children = crossover(parents[idx1], parents[idx2]);
            
            // Apply mutation
            Individual mutatedChild1 = mutate(children.first, MUTATION_RATE);
            Individual mutatedChild2 = mutate(children.second, MUTATION_RATE);
            
            // Add children to new population
            if (newPopulation.size() < population_.size()) {
                newPopulation.push_back(mutatedChild1);
            }
            if (newPopulation.size() < population_.size()) {
                newPopulation.push_back(mutatedChild2);
            }
        }
        
        // Replace old population with new population
        if (!newPopulation.empty()) {
            population_ = newPopulation;
        } else {
            std::cout << "Warning: New population is empty!" << std::endl;
        }
        
        // Evaluate fitness of the new population
        evaluateFitness();
    }
    
    // Find the best solution
    if (!population_.empty()) {
        bestSolution_ = findBestIndividual(population_);
        std::cout << "Genetic algorithm completed. Best solution found with fitness: " 
                  << bestSolution_.fitness << std::endl;
    } else {
        std::cout << "Error: Population is empty after evolution!" << std::endl;
    }
}

double ProcessManager::calculateFitness(Individual &individual) {
    // Initialize stocks and clear logs
    currentStocks_.clear();
    for (const auto &stock : config_.getStocks()) {
        currentStocks_[stock.name] = stock.quantity;
    }
    executionLogs_.clear();
    currentCycle_ = 0;

    // Track running processes with their completion times
    struct RunningProcess
    {
        const Process *process;
        int completionCycle;
    };
    vector<RunningProcess> runningProcesses;

    size_t processIndex = 0;
    while (currentCycle_ < delayLimit_ && processIndex < individual.processSequence.size())
    {
        // First, check if any processes have completed at this cycle
        auto it = runningProcesses.begin();
        while (it != runningProcesses.end())
        {
            if (it->completionCycle <= currentCycle_)
            {
                // Process completed, apply outputs
                for (const auto &output : it->process->outputs)
                {
                    currentStocks_[output.first] += output.second;
                }
                it = runningProcesses.erase(it);
            }
            else
            {
                ++it;
            }
        }

        // Try to execute the current process in the sequence
        string currentProcessName = individual.processSequence[processIndex];
        vector<const Process *> runnableProcesses = getRunnableProcesses();

        const Process *processToExecute = nullptr;
        for (const auto &proc : runnableProcesses)
        {
            if (proc->name == currentProcessName)
            {
                processToExecute = proc;
                break;
            }
        }

        if (processToExecute)
        {
            // Consume inputs immediately
            for (const auto &input : processToExecute->inputs)
            {
                currentStocks_[input.first] -= input.second;
            }

            // Schedule completion
            runningProcesses.push_back({processToExecute,
                                        currentCycle_ + processToExecute->nbCycle});

            // Log the start of the process
            executionLogs_.push_back({currentCycle_, processToExecute->name});

            // Move to next process in sequence
            processIndex++;
        }
        else
        {
            // Current process can't run, either due to resources or still running
            // Advance time to the next completion or just increment by 1
            if (!runningProcesses.empty())
            {
                int nextCompletionTime = delayLimit_;
                for (const auto &rp : runningProcesses)
                {
                    nextCompletionTime = min(nextCompletionTime, rp.completionCycle);
                }
                currentCycle_ = nextCompletionTime;
            }
            else
            {
                // No processes running, try next process
                processIndex++;
                currentCycle_++;
            }
        }
    }

    // Make sure all running processes complete (or reach the delay limit)
    while (!runningProcesses.empty() && currentCycle_ < delayLimit_)
    {
        int nextCompletionTime = delayLimit_;
        for (const auto &rp : runningProcesses)
        {
            nextCompletionTime = min(nextCompletionTime, rp.completionCycle);
        }
        currentCycle_ = nextCompletionTime;

        auto it = runningProcesses.begin();
        while (it != runningProcesses.end())
        {
            if (it->completionCycle <= currentCycle_)
            {
                // Process completed, apply outputs
                for (const auto &output : it->process->outputs)
                {
                    currentStocks_[output.first] += output.second;
                }
                it = runningProcesses.erase(it);
            }
            else
            {
                ++it;
            }
        }
    }

    // Calculate fitness as before
    double fitnessScore = 0.0;
    const vector<string>& optimizeGoals = config_.getOptimizeGoal();
    if (!optimizeGoals.empty()) {
        string goal = optimizeGoals[0];
        if (currentStocks_.count(goal) > 0)
            fitnessScore = currentStocks_[goal];
        else if (goal == "time")
            fitnessScore = -currentCycle_;
        else
            fitnessScore = numeric_limits<double>::lowest();
    } else {
        fitnessScore = 0.0;
    }
    return fitnessScore;
}

std::vector<Individual> ProcessManager::selectParents(const std::vector<Individual>& population) {
    std::vector<Individual> parents;
    if (population.empty()) {
        return parents;
    }
    
    // Find best individuals and total fitness
    double totalFitness = 0.0;
    double minFitness = std::numeric_limits<double>::max();
    
    for (const auto& individual : population) {
        totalFitness += individual.fitness;
        if (individual.fitness < minFitness) {
            minFitness = individual.fitness;
        }
    }
    
    // If all fitness values are negative or zero, shift them to positive range
    std::vector<double> adjustedFitness;
    if (totalFitness <= 0.0) {
        double shift = std::abs(minFitness) + 1.0;
        totalFitness = 0.0;
        for (const auto& individual : population) {
            double adjusted = individual.fitness + shift;
            adjustedFitness.push_back(adjusted);
            totalFitness += adjusted;
        }
    } else {
        for (const auto& individual : population) {
            adjustedFitness.push_back(individual.fitness);
        }
    }
    
    // Select parents using roulette wheel selection
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<> distribution(0.0, totalFitness);
    
    int numParentsToSelect = std::min(POPULATION_SIZE, static_cast<int>(population.size()));
    for (int i = 0; i < numParentsToSelect; ++i) {
        double randomValue = distribution(generator);
        double cumulativeProbability = 0.0;
        
        for (size_t j = 0; j < population.size(); ++j) {
            cumulativeProbability += adjustedFitness[j];
            if (randomValue <= cumulativeProbability) {
                parents.push_back(population[j]);
                break;
            }
        }
        
        // Ensure a parent is selected even if rounding errors occur
        if (parents.size() <= static_cast<size_t>(i) && !population.empty()) {
            parents.push_back(population[0]);
        }
    }
    
    return parents;
}

std::pair<Individual, Individual> ProcessManager::crossover(const Individual& parent1, const Individual& parent2) {
    if (parent1.processSequence.empty() || parent2.processSequence.empty()) {
        return std::make_pair(parent1, parent2); // Return copies of parents if either is empty
    }

    std::random_device rd;
    std::mt19937 generator(rd());
    
    // Calculate crossover point safely
    size_t minSize = std::min(parent1.processSequence.size(), parent2.processSequence.size());
    
    // If both sequences are very short, handle specially
    if (minSize == 0) {
        return std::make_pair(parent1, parent2);
    }
    
    // Generate crossover point between 0 and minSize-1
    std::uniform_int_distribution<> distribution(0, minSize - 1);
    size_t crossoverPoint = distribution(generator);
    
    // Create new sequences
    std::vector<std::string> child1Sequence;
    std::vector<std::string> child2Sequence;
    
    // First part of child1 comes from parent1
    for (size_t i = 0; i < crossoverPoint && i < parent1.processSequence.size(); ++i) {
        child1Sequence.push_back(parent1.processSequence[i]);
    }
    
    // Second part of child1 comes from parent2
    for (size_t i = crossoverPoint; i < parent2.processSequence.size(); ++i) {
        child1Sequence.push_back(parent2.processSequence[i]);
    }
    
    // First part of child2 comes from parent2
    for (size_t i = 0; i < crossoverPoint && i < parent2.processSequence.size(); ++i) {
        child2Sequence.push_back(parent2.processSequence[i]);
    }
    
    // Second part of child2 comes from parent1
    for (size_t i = crossoverPoint; i < parent1.processSequence.size(); ++i) {
        child2Sequence.push_back(parent1.processSequence[i]);
    }
    
    return std::make_pair(Individual(child1Sequence), Individual(child2Sequence));
}

vector<Individual> ProcessManager::crossoverPopulation(const vector<Individual>& parents, double crossoverRate) {
    vector<Individual> nextGeneration;
    random_device rd;
    mt19937 generator(rd());
    uniform_real_distribution<> probabilityDistribution(0.0, 1.0);
    for (int i = 0; i < POPULATION_SIZE / 2; ++i) {
        if (probabilityDistribution(generator) < crossoverRate) {
            pair<Individual, Individual> children = crossover(parents[i % parents.size()], parents[(i + 1) % parents.size()]);
            nextGeneration.push_back(children.first);
            nextGeneration.push_back(children.second);
        } else {
            nextGeneration.push_back(parents[i % parents.size()]);
            nextGeneration.push_back(parents[(i + 1) % parents.size()]);
        }
    }
    return nextGeneration;
}

Individual ProcessManager::mutate(const Individual& individual, double mutationRate) {
    // Get available processes
    std::vector<Process> availableProcesses = config_.getProcesses();
    if (availableProcesses.empty() || individual.processSequence.empty()) {
        return Individual(individual.processSequence);
    }
    
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<> probabilityDistribution(0.0, 1.0);
    
    std::vector<std::string> mutatedSequence = individual.processSequence;
    
    // For each position in the sequence, decide if it should mutate
    for (size_t i = 0; i < mutatedSequence.size(); ++i) {
        if (probabilityDistribution(generator) < mutationRate) {
            // Choose a random process
            std::uniform_int_distribution<> processDistribution(0, availableProcesses.size() - 1);
            int processIndex = processDistribution(generator);
            
            // Replace the process at position i
            mutatedSequence[i] = availableProcesses[processIndex].name;
        }
    }
    
    return Individual(mutatedSequence);
}

vector<Individual> ProcessManager::mutationPopulation(vector<Individual>& population, double mutationRate) {
    vector<Individual> mutatedPopulation;
    for (Individual& individual : population) {
        mutatedPopulation.push_back(mutate(individual, mutationRate));
    }
    return mutatedPopulation;
}

vector<const Process*> ProcessManager::getRunnableProcesses() {
    vector<const Process*> runnable;
    cout << "    getRunnableProcesses() - Current Stocks: ["; // Log current stocks at start
    for (const auto& pair : currentStocks_) {
        cout << pair.first << ":" << pair.second << ", ";
    }
    cout << "]" << endl;

    for (const auto &process : config_.getProcesses()) {
        cout << "    Checking process: " << process.name << endl; // Log process being checked
        bool canRun = true;
        for (const auto &input : process.inputs) {
            cout << "      Requires input: " << input.first << ":" << input.second << endl; // Log required input
            if (currentStocks_.count(input.first) == 0) {
                canRun = false;
                cout << "      Stock '" << input.first << "' not found in currentStocks_." << endl; // Log stock not found
                break;
            }
            if (currentStocks_[input.first] < input.second) {
                canRun = false;
                cout << "      Not enough stock '" << input.first << "'. Available: " << currentStocks_[input.first] << ", Required: " << input.second << endl; // Log insufficient stock
                break;
            }
        }
        if (canRun) {
            runnable.push_back(&process);
            cout << "    Process '" << process.name << "' is runnable." << endl; // Log runnable process
        } else {
            cout << "    Process '" << process.name << "' is NOT runnable." << endl; // Log non-runnable process
        }
    }
    cout << "    getRunnableProcesses() - Runnable processes count: " << runnable.size() << endl; // Log count of runnable processes
    return runnable;
}

bool ProcessManager::executeProcess(const Process* process) {
    if (!process)
        return false;
    for (const auto &input : process->inputs) {
        currentStocks_[input.first] -= input.second;
        if (currentStocks_[input.first] < 0)
            return false;
    }
    for (const auto &output : process->outputs) {
        currentStocks_[output.first] += output.second;
    }
    return true;
}

void ProcessManager::generateOutput() {
    // Find the best individual
    double bestFitness = numeric_limits<double>::lowest();
    const Individual *bestIndividual = nullptr;
    for (const auto &indiv : population_)
    {
        if (indiv.fitness > bestFitness)
        {
            bestFitness = indiv.fitness;
            bestIndividual = &indiv;
        }
    }

    if (!bestIndividual)
    {
        cout << "\nNo best individual found in final population." << endl;
        return;
    }

    // Re-run the simulation with the best individual to generate accurate logs
    cout << "\nRe-running simulation with best individual..." << endl;

    // Clear stocks and logs
    currentStocks_.clear();
    for (const auto &stock : config_.getStocks())
    {
        currentStocks_[stock.name] = stock.quantity;
    }
    executionLogs_.clear();
    currentCycle_ = 0;

    // Simulate the best individual's process sequence
    struct RunningProcess
    {
        const Process *process;
        int completionCycle;
    };
    vector<RunningProcess> runningProcesses;

    size_t processIndex = 0;
    while (currentCycle_ < delayLimit_ && processIndex < bestIndividual->processSequence.size())
    {
        // Check for completed processes at this cycle
        auto it = runningProcesses.begin();
        while (it != runningProcesses.end())
        {
            if (it->completionCycle <= currentCycle_)
            {
                // Process completed, apply outputs
                for (const auto &output : it->process->outputs)
                {
                    currentStocks_[output.first] += output.second;
                }
                it = runningProcesses.erase(it);
            }
            else
            {
                ++it;
            }
        }

        // Try to execute the current process in the sequence
        const string &currentProcessName = bestIndividual->processSequence[processIndex];
        vector<const Process *> runnableProcesses = getRunnableProcesses();

        const Process *processToExecute = nullptr;
        for (const auto &proc : runnableProcesses)
        {
            if (proc->name == currentProcessName)
            {
                processToExecute = proc;
                break;
            }
        }

        if (processToExecute)
        {
            // Consume inputs
            for (const auto &input : processToExecute->inputs)
            {
                currentStocks_[input.first] -= input.second;
            }

            // Schedule completion
            runningProcesses.push_back({processToExecute,
                                        currentCycle_ + processToExecute->nbCycle});

            // Log the process start
            executionLogs_.push_back({currentCycle_, processToExecute->name});

            processIndex++;
        }
        else
        {
            // Process can't run yet or ever
            if (!runningProcesses.empty())
            {
                // Advance to next completion
                int nextCompletionTime = delayLimit_;
                for (const auto &rp : runningProcesses)
                {
                    nextCompletionTime = min(nextCompletionTime, rp.completionCycle);
                }
                currentCycle_ = nextCompletionTime;
            }
            else
            {
                // Skip to next process or advance time
                processIndex++;
                currentCycle_++;
            }
        }
    }

    // Complete any remaining running processes
    while (!runningProcesses.empty() && currentCycle_ < delayLimit_)
    {
        int nextCompletionTime = delayLimit_;
        for (const auto &rp : runningProcesses)
        {
            nextCompletionTime = min(nextCompletionTime, rp.completionCycle);
        }
        currentCycle_ = nextCompletionTime;

        auto it = runningProcesses.begin();
        while (it != runningProcesses.end())
        {
            if (it->completionCycle <= currentCycle_)
            {
                for (const auto &output : it->process->outputs)
                {
                    currentStocks_[output.first] += output.second;
                }
                it = runningProcesses.erase(it);
            }
            else
            {
                ++it;
            }
        }
    }

    // Now display accurate results in the required format
    cout << "Nice file! " << config_.getProcesses().size() << " processes, "
         << config_.getStocks().size() << " stocks, "
         << config_.getOptimizeGoal().size() << " to optimize" << endl;

    cout << "Evaluating .................. done." << endl;
    cout << "Main walk" << endl;

    // Print process execution log in the required format
    for (auto log : executionLogs_)
    {
        cout << log.first << ":" << log.second << endl;
    }

    cout << "no more process doable at time " << currentCycle_ << endl;

    cout << "Stock :" << endl;
    for (auto stock : currentStocks_)
    {
        cout << stock.first << "=> " << stock.second << endl;
    }
}
