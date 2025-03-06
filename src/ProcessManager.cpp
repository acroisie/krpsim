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
        std::vector<const Process*> runnableProcesses = getRunnableProcesses();
        if (runnableProcesses.empty()) {
            std::cout << "No more processes to run." << std::endl;
            break;
        }

        const Process* choosenProcess = chooseGreedyProcesses(runnableProcesses);
        if (choosenProcess) {
            executeProcess(choosenProcess);
        } else {
            std::cout << "No process to execute." << std::endl;
            break;
        }

        currentCycle_++;
        updateStocksWithOutputs();
    }
    generateOutput();
    return true;
}

std::vector<const Process*> ProcessManager::getRunnableProcesses() {
    std::vector<const Process*> runnable;
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

const Process* ProcessManager::chooseGreedyProcesses(const std::vector<const Process*>& runnableProcesses) {
    if (runnableProcesses.empty()) {
        return nullptr;
    }
	
    return runnableProcesses[0];
}

bool ProcessManager::executeProcess(const Process* process) {
    if (!process) {
        return false;
    }

    for (const auto &input : process->inputs) {
        currentStocks_[input.first] -= input.second;
        if (currentStocks_[input.first] < 0) {
            std::cout << "Error: Not enough input stock " << input.first 
                      << " for process " << process->name << std::endl;
            return false;
        }
    }

    for (const auto &output : process->outputs) {
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