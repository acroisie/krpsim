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

const int ProcessManager::POPULATION_SIZE = 50;

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

bool ProcessManager::runGeneticAlgorithm() {
    double mutationRate = 0.1;
    double crossoverRate = 0.7;
    int generationCount = 50;
    cout << "Running genetic algorithm... (runGeneticAlgorithm() START)" << endl;
    initializePopulation();
    cout << "Initial population generated." << endl;
    cout << "Starting evolution over " << generationCount << " generations..." << endl;

    // Evaluate initial population fitness
    evaluateFitness();

    for (int generation = 0; generation < generationCount; ++generation) {
        cout << "--- Generation " << generation << " ---" << endl;

        // Selection, crossover, and mutation
        vector<Individual> parents = selection();
        vector<Individual> nextGeneration = crossoverPopulation(parents, crossoverRate);
        population_ = mutationPopulation(nextGeneration, mutationRate);

        // IMPORTANT FIX: Evaluate the new population's fitness
        evaluateFitness();

        // Now find the best fitness in the properly evaluated population
        double bestFitnessThisGen = numeric_limits<double>::lowest();
        for (const auto& indiv : population_) {
            if (indiv.fitness > bestFitnessThisGen)
                bestFitnessThisGen = indiv.fitness;
        }
        cout << "  Best fitness in generation " << generation << ": " << bestFitnessThisGen << endl;
    }

    cout << "Evolution finished over " << generationCount << " generations." << endl;
    generateOutput();
    cout << "Running genetic algorithm... (runGeneticAlgorithm() END)" << endl;
    return true;
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

vector<Individual> ProcessManager::selectParents(const vector<Individual>& population) {
    vector<Individual> parents;
    vector<double> fitnessValues;
    for (const auto& individual : population) {
        fitnessValues.push_back(individual.fitness);
    }
    double totalFitness = 0.0;
    for (double fitness : fitnessValues) {
        totalFitness += fitness;
    }
    vector<double> selectionProbabilities;
    if (totalFitness <= 0)
        selectionProbabilities.assign(fitnessValues.size(), 1.0 / fitnessValues.size());
    else
    {
        for (double fitness : fitnessValues)
        {
            selectionProbabilities.push_back(fitness / totalFitness);
        }
    }
    random_device rd;
    mt19937 generator(rd());
    uniform_real_distribution<> distribution(0.0, 1.0);
    for (int i = 0; i < POPULATION_SIZE; ++i) {
        double randomValue = distribution(generator);
        double cumulativeProbability = 0.0;
        for (size_t j = 0; j < selectionProbabilities.size(); ++j) {
            cumulativeProbability += selectionProbabilities[j];
            if (randomValue <= cumulativeProbability) {
                parents.push_back(population[j]);
                break;
            }
        }
    }
    return parents;
}

pair<Individual, Individual> ProcessManager::crossover(const Individual& parent1, const Individual& parent2) {
    random_device rd;
    mt19937 generator(rd());
    uniform_int_distribution<> distribution(0, min(parent1.processSequence.size(), parent2.processSequence.size()) > 0 ? min(parent1.processSequence.size(), parent2.processSequence.size()) - 1 : 0);
    size_t crossoverPoint = 0;
    if (parent1.processSequence.size() > 0 && parent2.processSequence.size() > 0)
        crossoverPoint = distribution(generator);
    vector<string> child1Sequence;
    vector<string> child2Sequence;
    for (size_t i = 0; i < parent1.processSequence.size(); ++i) {
        if (i <= crossoverPoint)
            child1Sequence.push_back(parent1.processSequence[i]);
        else if (i < parent2.processSequence.size())
            child1Sequence.push_back(parent2.processSequence[i]);
    }
    for (size_t i = 0; i < parent2.processSequence.size(); ++i) {
        if (i <= crossoverPoint && i < parent1.processSequence.size())
            child2Sequence.push_back(parent2.processSequence[i]);
         else if (i < parent1.processSequence.size())
             child2Sequence.push_back(parent1.processSequence[i]);
    }
    return make_pair(Individual(child1Sequence), Individual(child2Sequence));
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
    random_device rd;
    mt19937 generator(rd());
    uniform_real_distribution<> probabilityDistribution(0.0, 1.0);
    uniform_int_distribution<> processDistribution(0, config_.getProcesses().size() - 1);
    vector<string> mutatedSequence = individual.processSequence;
    if (probabilityDistribution(generator) < mutationRate) {
        if (!mutatedSequence.empty()) {
            uniform_int_distribution<> pointDistribution(0, mutatedSequence.size() - 1);
            int mutationPoint = pointDistribution(generator);
            int processIndex = processDistribution(generator);
            string processNameToMutateTo = config_.getProcesses()[processIndex].name;
            mutatedSequence[mutationPoint] = processNameToMutateTo;
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
    for (const auto &process : config_.getProcesses()) {
        bool canRun = true;
        for (const auto &input : process.inputs) {
            if (currentStocks_.count(input.first) == 0) {
                canRun = false;
                break;
            }
            if (currentStocks_[input.first] < input.second) {
                canRun = false;
                break;
            }
        }
        if (canRun)
            runnable.push_back(&process);
    }
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
