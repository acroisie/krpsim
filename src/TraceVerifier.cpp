#include "TraceVerifier.hpp"
#include "StringUtils.hpp"
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

TraceVerifier::TraceVerifier(const Config &configRef, int timeLimit)
    : config(configRef), timeLimit(timeLimit) {
    for (const auto &stock : config.getStocks())
        stocks[stock.name] = stock.quantity;
}

bool TraceVerifier::canStart(const Process *process) const {
    if (!process) return false;
    for (const auto &[resource, quantity] : process->inputs) {
        auto stockIt = stocks.find(resource);
        if (stockIt == stocks.end() || stockIt->second < quantity) return false;
    }
    return true;
}

void TraceVerifier::addOutputs(const Process *process) {
    if (!process) return;
    for (const auto &[resource, quantity] : process->outputs)
        stocks[resource] += quantity;
}

void TraceVerifier::finishUntil(int cycle) {
    while (!runningQ.empty() && runningQ.top().completionTime <= cycle) {
        addOutputs(runningQ.top().proc);
        runningQ.pop();
    }
}

bool TraceVerifier::start(int cycle, const Process *process, string &error) {
    if (!process) {
        error = "process inconnu";
        return false;
    }
    if (cycle > timeLimit) {
        error = "dépasse la limite de temps";
        return false;
    }
    if (!canStart(process)) {
        error = "stock insuffisant";
        return false;
    }

    for (const auto &[resource, quantity] : process->inputs)
        stocks[resource] -= quantity;
    runningQ.push({process, cycle + process->nbCycle});
    return true;
}

void TraceVerifier::printStocks(int finalCycle) const {
    cout << "Stock :" << endl;
    set<string> stockNames;
    for (const auto &stock : config.getStocks()) stockNames.insert(stock.name);
    for (const auto &process : config.getProcesses()) {
        for (const auto &[resource, _] : process.inputs)
            stockNames.insert(resource);
        for (const auto &[resource, _] : process.outputs)
            stockNames.insert(resource);
    }
    for (const string &name : stockNames)
        cout << "  " << name << " => "
             << (stocks.count(name) ? stocks.at(name) : 0) << endl;
    cout << "----------------------------------------" << endl;
    cout << "Dernier cycle : " << finalCycle << endl;
}

bool TraceVerifier::verifyFile(const string &traceFile) {
    ifstream file(traceFile);
    if (!file) {
        cerr << "Erreur : impossible d’ouvrir " << traceFile << endl;
        return false;
    }

    string line;
    int previousCycle = 0, lastCycle = 0;
    size_t lineNumber = 0;

    while (getline(file, line)) {
        ++lineNumber;
        line = StringUtils::trim(line);
        if (line.empty() || line[0] == '#') continue;

        size_t colonIndex = line.find(':');
        if (colonIndex == string::npos) {
            cerr << "Syntaxe invalide ligne " << lineNumber << endl;
            return false;
        }

        int cycle = stoi(line.substr(0, colonIndex));
        string processName = StringUtils::trim(line.substr(colonIndex + 1));
        if (cycle < previousCycle) {
            cerr << "Cycles non croissants ligne " << lineNumber << endl;
            return false;
        }

        finishUntil(cycle);

        const Process *processPtr = nullptr;
        for (const auto &process : config.getProcesses())
            if (process.name == processName) {
                processPtr = &process;
                break;
            }

        string error;
        if (!start(cycle, processPtr, error)) {
            cerr << "Erreur ligne " << lineNumber << " (" << cycle << ":"
                 << processName << ") ➜ " << error << endl;
            return false;
        }

        previousCycle = cycle;
        lastCycle = cycle;
    }

    if (!runningQ.empty())
        lastCycle = max(lastCycle, runningQ.top().completionTime);
    finishUntil(lastCycle);

    cout << "Trace vérifiée avec succès." << endl
         << "----------------------------------------" << endl;
    printStocks(lastCycle);
    return true;
}
