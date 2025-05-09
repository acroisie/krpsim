#include "TraceVerifier.hpp"
#include "StringUtils.hpp"
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

namespace {
    void addProcessOutputs(const Process *process, map<string, int> &stocks) {
        if (!process) return;
        for (const auto &[resource, quantity] : process->outputs)
            stocks[resource] += quantity;
    }

    const Process *findProcess(const vector<Process> &procs,
                               const string &name) {
        for (const auto &p : procs)
            if (p.name == name) return &p;
        return nullptr;
    }
} // namespace

TraceVerifier::TraceVerifier(const Config &cfg, int limit)
    : config(cfg), timeLimit(limit) {
    for (const auto &stock : config.getStocks())
        stocks[stock.name] = stock.quantity;
}

bool TraceVerifier::canStart(const Process *process) const {
    if (!process) return false;
    for (const auto &[resource, quantity] : process->inputs) {
        auto it = stocks.find(resource);
        if (it == stocks.end() || it->second < quantity) return false;
    }
    return true;
}

void TraceVerifier::finishUntil(int cycle) {
    while (!runningQ.empty() && runningQ.top().completionTime <= cycle) {
        addProcessOutputs(runningQ.top().proc, stocks);
        runningQ.pop();
    }
}

bool TraceVerifier::start(int cycle, const Process *process, string &err) {
    if (!process) {
        err = "process inconnu";
        return false;
    }
    if (cycle > timeLimit) {
        err = "dépasse la limite de temps";
        return false;
    }
    if (!canStart(process)) {
        err = "stock insuffisant";
        return false;
    }
    for (const auto &[resource, quantity] : process->inputs)
        stocks[resource] -= quantity;
    runningQ.push({process, cycle + process->nbCycle});
    return true;
}

void TraceVerifier::printStocks(int finalCycle) const {
    cout << "Stock :" << endl;
    set<string> names;
    for (const auto &stock : config.getStocks()) names.insert(stock.name);
    for (const auto &process : config.getProcesses()) {
        for (const auto &[resource, _] : process.inputs) names.insert(resource);
        for (const auto &[resource, _] : process.outputs)
            names.insert(resource);
    }
    for (const string &name : names)
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
    int prev = 0, last = 0;
    size_t nline = 0;
    while (getline(file, line)) {
        ++nline;
        line = StringUtils::trim(line);
        if (line.empty() || line[0] == '#') continue;
        size_t col = line.find(':');
        if (col == string::npos) {
            cerr << "Syntaxe invalide ligne " << nline << endl;
            return false;
        }
        int cycle = stoi(line.substr(0, col));
        string name = StringUtils::trim(line.substr(col + 1));
        if (cycle < prev) {
            cerr << "Cycles non croissants ligne " << nline << endl;
            return false;
        }
        finishUntil(cycle);
        const Process *p = findProcess(config.getProcesses(), name);
        string err;
        if (!start(cycle, p, err)) {
            cerr << "Erreur ligne " << nline << " (" << cycle << ":" << name
                 << ") ➜ " << err << endl;
            return false;
        }
        prev = cycle;
        last = cycle;
    }
    if (!runningQ.empty()) last = max(last, runningQ.top().completionTime);
    finishUntil(last);
    cout << "Trace vérifiée avec succès." << endl
         << "----------------------------------------" << endl;
    printStocks(last);
    return true;
}
