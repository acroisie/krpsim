#include "ProcessManager.hpp"
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

bool ProcessManager::runSimulation() {
    while (currentCycle_ < delayLimit_) {
        vector<const Process*> runnableProcesses = getRunnableProcesses();
        if (runnableProcesses.empty()) {
            cout << "No more processes to run." << endl;
            break;
        }
        
        currentCycle_++;
    }
    generateOutput();
    return true;
}

bool ProcessManager::runGeneticAlgorithm() {
    cout << "Running genetic algorithm... (runGeneticAlgorithm() START)" << endl; // DEBUG: runGeneticAlgorithm start

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

        cout << "Individual " << i << ": Sequence = [ ";
        for (const auto& processName : population_[i].processSequence) {
            cout << processName << ", ";
        }
        cout << "]" << endl;
    }

    for (int i = 0; i < POPULATION_SIZE; ++i) {
        population_[i].fitness = calculateFitness(population_[i]);
        cout << "Individual " << i << ": Fitness = " << population_[i].fitness << endl;
    }

    generateOutput();
    cout << "Running genetic algorithm... (runGeneticAlgorithm() END)" << endl; // DEBUG: runGeneticAlgorithm end
    return true;
}

double ProcessManager::calculateFitness(Individual &individual) {
    cout << "calculateFitness() called for individual. Sequence: [";
    for (const auto& procName : individual.processSequence) cout << procName << ", ";
    cout << "]" << endl;

    currentStocks_.clear();
    for (const auto &stock : config_.getStocks()) {
        currentStocks_[stock.name] = stock.quantity;
    }
    cout << "  Initial stocks in calculateFitness: ["; // Ajout du log pour l'initialisation
    for (const auto& stockItem : currentStocks_) cout << stockItem.first << ":" << stockItem.second << ", ";
    cout << "]" << endl;

    executionLogs_.clear();
    currentCycle_ = 0;

    for (const string& processName : individual.processSequence) {
        if (currentCycle_ >= delayLimit_) {
            break;
        }

        cout << "  Cycle " << currentCycle_ << ": Starting cycle, process from sequence: " << processName << endl;

        bool processExecutedThisCycle = false;

        do {
            processExecutedThisCycle = false;
            cout << "    Cycle " << currentCycle_ << ": Entering do-while loop." << endl;
            vector<const Process*> runnableProcesses = getRunnableProcesses();
            cout << "    Cycle " << currentCycle_ << ": getRunnableProcesses() returned " << runnableProcesses.size() << " processes." << endl;

            if (!runnableProcesses.empty()) {
                const Process* processToExecute = nullptr;
                for (const auto &proc : runnableProcesses) {
                    if (proc->name == processName) {
                        processToExecute = proc;
                        break;
                    }
                }
                if (!processToExecute && !runnableProcesses.empty()) {
                    processToExecute = runnableProcesses[0];
                }


                if (processToExecute) {
                    cout << "    Cycle " << currentCycle_ << ": Executing process: " << processToExecute->name << endl;
                    executeProcess(processToExecute);
                    executionLogs_.push_back({currentCycle_, processToExecute->name});
                    processExecutedThisCycle = true;
                } else {
                    cout << "    Cycle " << currentCycle_ << ": Aucun processus de sequence executable et aucun autre processus runnable trouvé." << endl;
                    break;
                }
            } else {
                 cout << "    Cycle " << currentCycle_ << ": Aucun processus exécutable trouvé par getRunnableProcesses()." << endl;
                break;
            }

        } while (processExecutedThisCycle);

        currentCycle_++;
    }

    double fitnessScore = 0.0;
    const vector<string>& optimizeGoals = config_.getOptimizeGoal();
    if (!optimizeGoals.empty()) {
        string goal = optimizeGoals[0];
        if (currentStocks_.count(goal) > 0) {
            fitnessScore = currentStocks_[goal];
        } else if (goal == "time") {
            fitnessScore = -currentCycle_;
        } else {
            cout << "Warning: Optimize goal '" << goal << "' not found in stocks." << endl;
            fitnessScore = numeric_limits<double>::lowest();
        }
    } else {
        fitnessScore = 0.0;
    }
    cout << "calculateFitness() finished, fitnessScore = " << fitnessScore << endl; // DEBUG: Fin calculateFitness et fitness
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
    for (double fitness : fitnessValues) {
        selectionProbabilities.push_back(fitness / totalFitness);
    }

    random_device rd;
    mt19937 generator(rd());
    uniform_real_distribution<> distribution(0.0, 1.0);

    for (int i = 0; i < 2; ++i) {
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
    uniform_int_distribution<> distribution(0, parent1.processSequence.size() - 1);

    size_t crossoverPoint = distribution(generator);
    vector<string> child1Sequence;
    vector<string> child2Sequence;

    for (size_t i = 0; i < parent1.processSequence.size(); ++i) {
        if (i <= crossoverPoint) {
            child1Sequence.push_back(parent1.processSequence[i]);
            child2Sequence.push_back(parent2.processSequence[i]);
        } else { 
            child1Sequence.push_back(parent2.processSequence[i]);
            child2Sequence.push_back(parent1.processSequence[i]);
        }
    }

    return make_pair(Individual(child1Sequence), Individual(child2Sequence));
}

Individual ProcessManager::mutate(const Individual& individual) {
    random_device rd;
    mt19937 generator(rd());
    uniform_int_distribution<> processDistribution(0, individual.processSequence.size() - 1);
    uniform_int_distribution<> lengthDistribution(1, 3);

    vector<string> mutatedSequence = individual.processSequence;
    size_t mutationPoint = processDistribution(generator);
    size_t mutationLength = lengthDistribution(generator);

    for (size_t i = mutationPoint; i < mutationPoint + mutationLength; ++i) {
        if (i < mutatedSequence.size()) {
            mutatedSequence[i] = "mutated"; // Adapt the rest of the code 
        }
    }

    return Individual(mutatedSequence);
}

vector<const Process*> ProcessManager::getRunnableProcesses() {
    cout << "  getRunnableProcesses() called, cycle: " << currentCycle_ << ", currentStocks: [";
    for (const auto& stockItem : currentStocks_) cout << stockItem.first << ":" << stockItem.second << ", ";
    cout << "]" << endl;

    vector<const Process*> runnable;
    for (const auto &process : config_.getProcesses()) {
        bool canRun = true;
        cout << "    Checking process: " << process.name << ", inputs: [";
        for (const auto &input : process.inputs) {
            cout << input.first << ":" << input.second << ", ";
            if (currentStocks_.count(input.first) == 0) { // Vérifier si le stock existe dans currentStocks_
                canRun = false;
                cout << "] - NOT RUNNABLE, stock manquant: " << input.first << " (requis mais non défini dans currentStocks_)" << endl;
                break;
            }
            if (currentStocks_[input.first] < input.second) {
                canRun = false;
                cout << "] - NOT RUNNABLE, stock insuffisant: " << input.first << " (requis: " << input.second << ", disponible: " << currentStocks_[input.first] << ")" << endl;
                break;
            }
        }
        if (canRun) {
            cout << "] - RUNNABLE" << endl;
            runnable.push_back(&process);
        }
    }
    cout << "  getRunnableProcesses() returning " << runnable.size() << " processes." << endl;
    return runnable;
}

bool ProcessManager::executeProcess(const Process* process) {
    if (!process) {
        return false;
    }

    cout << "Executing process: " << process->name << endl;

    for (const auto &input : process->inputs) {
        cout << "  Consuming stock: '" << input.first << "', quantity: " << input.second << endl;
        currentStocks_[input.first] -= input.second;
        if (currentStocks_[input.first] < 0) {
            cout << "Error: Not enough input stock " << input.first
                      << " for process " << process->name << endl;
            return false;
        }
    }

    for (const auto &output : process->outputs) {
        cout << "  Producing stock: '" << output.first << "', quantity: " << output.second << endl;
        currentStocks_[output.first] += output.second;
    }

    return true;
}

void ProcessManager::generateOutput() {
	cout << "Simulation completed in " << currentCycle_ << " cycles." << endl;
	cout << "Final stocks:" << endl;
	for (auto stock : currentStocks_) {
		cout << "Stock " << stock.first << ": " << stock.second << endl;
	}
	cout << "Execution logs:" << endl;
	for (auto log : executionLogs_) {
		cout << "Cycle " << log.first << ": " << log.second << endl;
	}
}
