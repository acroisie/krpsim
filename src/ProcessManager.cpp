#include "ProcessManager.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <set>

using namespace std;

ProcessManager::ProcessManager(const Config &config, int timeLimit)
    : config(config), simulator(config, timeLimit),
      geneticAlgorithm(config, simulator, 100, 0.05, 0.7, 4, 50, 200),
      timeLimit(timeLimit) {}

void ProcessManager::run() {
    cout << "Nice file! " << config.getProcesses().size() << " processes, "
         << config.getStocks().size() << " initial stocks, "
         << config.getOptimizeGoal().size() << " optimization goal(s)" << endl;

    cout << "Evaluating using Genetic Algorithm..." << endl;

    int baseGenerations = 100;
    int processCount = config.getProcesses().size();
    int stockCount = config.getStocks().size();

    int generations = baseGenerations;
    if (processCount > 10 || stockCount > 5) {
        generations = 150;
    }
    if (processCount > 15 || stockCount > 10) {
        generations = 200;
    }

    Individual bestSolution = geneticAlgorithm.runEvolution(generations);

    cout << "Optimization complete." << endl;
    generateOutput(bestSolution);
}

void ProcessManager::generateOutput(const Individual &bestSolution) {
    cout << "----------------------------------------" << endl;
    cout << "Best solution found:" << endl;

    Simulator::SimulationResult result =
        simulator.runSimulation(bestSolution.processSequence);

    cout << "Final Fitness Score: " << result.fitness << endl;

    if (result.executionLog.empty()) {
        cout << "(No processes executed)" << endl;
    } else {
        for (const auto &[cycle, processName] : result.executionLog) {
            cout << cycle << ":" << processName << endl;
        }
    }

    const char *tracefile_env = std::getenv("KRPSIM_TRACEFILE");
    std::string tracefile = tracefile_env ? tracefile_env : "tracefile.txt";
    std::ofstream trace(tracefile);
    if (trace.is_open()) {
        for (const auto &[cycle, processName] : result.executionLog) {
            trace << cycle << ":" << processName << "\n";
        }
        trace.close();
        cout << "(Tracefile generated: " << tracefile << ")" << endl;
    } else {
        cerr << "Warning: Could not write tracefile '" << tracefile << "'"
             << endl;
    }

    cout << "----------------------------------------" << endl;

    if (result.executionLog.empty()) {
        cout << "No process could be executed within the time limit ("
             << timeLimit << ")." << endl;
    } else if (!result.timeoutReached && result.finalCycle < timeLimit) {
        cout << "no more process doable at time " << result.finalCycle << endl;
    } else {
        cout << "Simulation reached time limit at cycle " << timeLimit << "."
             << endl;
    }

    cout << "Stock:" << endl;

    set<string> relevantStocks;

    for (const auto &stock : config.getStocks()) {
        relevantStocks.insert(stock.name);
    }

    for (const auto &process : config.getProcesses()) {
        for (const auto &[resource, _] : process.inputs) {
            relevantStocks.insert(resource);
        }
        for (const auto &[resource, _] : process.outputs) {
            relevantStocks.insert(resource);
        }
    }

    for (const string &stockName : relevantStocks) {
        int quantity = 0;
        auto it = result.finalStocks.find(stockName);
        if (it != result.finalStocks.end()) {
            quantity = it->second;
        }
        cout << "  " << stockName << " => " << quantity << endl;
    }

    cout << "----------------------------------------" << endl;
}