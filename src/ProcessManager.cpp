#include "ProcessManager.hpp"
#include <iostream>
#include <limits>

using namespace std;

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

        // const Process* choosenProcess = chooseGreedyProcesses(runnableProcesses);
        // if (choosenProcess) {
        //     executeProcess(choosenProcess);
        // } else {
        //     cout << "No process to execute." << endl;
        //     break;
        // }

        currentCycle_++;
        updateStocksWithOutputs();
    }
    generateOutput();
    return true;
}

bool ProcessManager::runGeneticAlgorithm() {
    cout << "Genetic algorithm is not implemented yet." << endl;
    return true;
}

vector<const Process*> ProcessManager::getRunnableProcesses() {
    vector<const Process*> runnable;
    for (const auto &process : config_.getProcesses()) {
        bool canRun = true;
        for (const auto &input : process.inputs) {
            if (currentStocks_[input.first] < input.second) {
                canRun = false;
                break;
            }
        }
        if (canRun) {
            runnable.push_back(&process);
        }
    }
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

void ProcessManager::updateStocksWithOutputs() {
	// Nothing to do here
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