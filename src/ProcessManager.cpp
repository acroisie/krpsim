#include "ProcessManager.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <set>
using namespace std;

ProcessManager::ProcessManager(const Config &config, int timeLimit)
    : config(config), simulator(config, timeLimit),
      geneticAlgorithm(config, simulator, 100, 0.05, 0.7, 4, 50,
                       calculateMaxSequenceLength(config, timeLimit)),
      timeLimit(timeLimit) {}

int ProcessManager::calculateMaxSequenceLength(const Config &config,
                                               int timeLimit) {
    const auto &processes = config.getProcesses();
    const auto &stocks = config.getStocks();

    if (processes.empty()) {
        return 100;
    }

    double avgCycle = 0.0;
    for (const auto &process : processes) {
        avgCycle += process.nbCycle;
    }
    avgCycle /= processes.size();

    int estimatedProcesses =
        static_cast<int>(timeLimit / std::max(1.0, avgCycle));

    double complexityFactor =
        2.0 + (stocks.size() * 0.1) + (processes.size() * 0.2);

    int maxLength = static_cast<int>(estimatedProcesses * complexityFactor);

    int minValue = std::max(100, estimatedProcesses);
    int maxValue = std::min(20000, estimatedProcesses * 10);

    return std::max(minValue, std::min(maxValue, maxLength));
}

void ProcessManager::run() {
    cout << "Nice file! " << config.getProcesses().size() << " processes, "
         << config.getStocks().size() << " initial stocks, "
         << config.getOptimizeGoal().size() << " optimization goal(s)" << endl;
    cout << "Evaluating using Genetic Algorithm..." << endl;
    int processCount = config.getProcesses().size();
    int stockCount = config.getStocks().size();
    int generations = 100 + (processCount > 10 || stockCount > 5 ? 50 : 0) +
                      (processCount > 15 || stockCount > 10 ? 50 : 0);
    Individual bestIndividual = geneticAlgorithm.runEvolution(generations);
    cout << "Optimization complete." << endl;
    generateOutput(bestIndividual);
}

void ProcessManager::generateOutput(const Individual &bestIndividual) {
    cout << "----------------------------------------\nBest solution found:"
         << endl;
    auto result = simulator.runSimulation(bestIndividual.processSequence);
    cout << "Final Fitness Score: " << result.fitness << endl;
    if (result.executionLog.empty())
        cout << "(No processes executed)" << endl;
    else
        for (const auto &[cycle, processName] : result.executionLog)
            cout << cycle << ":" << processName << endl;
    const char *tracefileEnv = getenv("KRPSIM_TRACEFILE");
    string tracefile = tracefileEnv ? tracefileEnv : "tracefile.txt";
    ofstream trace(tracefile);
    if (trace.is_open()) {
        for (const auto &[cycle, processName] : result.executionLog)
            trace << cycle << ":" << processName << "\n";
        cout << "(Tracefile generated: " << tracefile << ")" << endl;
    } else
        cerr << "Warning: Could not write tracefile '" << tracefile << "'"
             << endl;
    cout << "----------------------------------------" << endl;
    if (result.executionLog.empty())
        cout << "No process could be executed within the time limit ("
             << timeLimit << ")." << endl;
    else if (!result.timeoutReached && result.finalCycle < timeLimit)
        cout << "No more process doable at time " << result.finalCycle << endl;
    else
        cout << "Simulation reached time limit at cycle " << timeLimit << "."
             << endl;
    cout << "Stock summary:" << endl;
    set<string> allStockNames;
    for (const auto &stock : config.getStocks())
        allStockNames.insert(stock.name);
    for (const auto &process : config.getProcesses()) {
        for (const auto &[resource, _] : process.inputs)
            allStockNames.insert(resource);
        for (const auto &[resource, _] : process.outputs)
            allStockNames.insert(resource);
    }
    for (const string &stockName : allStockNames)
        cout << "  " << stockName << " => "
             << (result.finalStocks.count(stockName)
                     ? result.finalStocks.at(stockName)
                     : 0)
             << endl;
    cout << "----------------------------------------" << endl;
}