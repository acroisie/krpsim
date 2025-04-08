#include "TraceVerifier.hpp"
#include "StringUtils.hpp"
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

TraceVerifier::TraceVerifier(const Config &cfg, int limit)
    : config(cfg), timeLimit(limit) {
    for (const auto &s : config.getStocks()) stocks[s.name] = s.quantity;
}

bool TraceVerifier::canStart(const Process *p) const {
    if (!p) return false;
    for (const auto &[r, q] : p->inputs) {
        auto it = stocks.find(r);
        if (it == stocks.end() || it->second < q) return false;
    }
    return true;
}

void TraceVerifier::addOutputs(const Process *p) {
    if (!p) return;
    for (const auto &[r, q] : p->outputs) stocks[r] += q;
}

void TraceVerifier::finishUntil(int cycle) {
    while (!runningQ.empty() && runningQ.top().completionTime <= cycle) {
        addOutputs(runningQ.top().proc);
        runningQ.pop();
    }
}

bool TraceVerifier::start(int cycle, const Process *p, string &err) {
    if (!p)         { err = "process inconnu";      return false; }
    if (cycle > timeLimit) { err = "dépasse la limite de temps"; return false; }
    if (!canStart(p)){ err = "stock insuffisant";   return false; }

    for (const auto &[r, q] : p->inputs) stocks[r] -= q;
    runningQ.push({p, cycle + p->nbCycle});
    return true;
}

void TraceVerifier::printStocks(int finalCycle) const {
    cout << "Stock :" << endl;
    set<string> names;
    for (const auto &s : config.getStocks()) names.insert(s.name);
    for (const auto &p : config.getProcesses()) {
        for (const auto &[r, _] : p.inputs)  names.insert(r);
        for (const auto &[r, _] : p.outputs) names.insert(r);
    }
    for (const string &n : names)
        cout << "  " << n << " => "
             << (stocks.count(n) ? stocks.at(n) : 0) << endl;
    cout << "----------------------------------------" << endl;
    cout << "Dernier cycle : " << finalCycle << endl;
}

bool TraceVerifier::verifyFile(const string &traceFile) {
    ifstream file(traceFile);
    if (!file) { cerr << "Erreur : impossible d’ouvrir " << traceFile << endl; return false; }

    string line;
    int prev = 0, last = 0;
    size_t nline = 0;

    while (getline(file, line)) {
        ++nline;
        line = StringUtils::trim(line);
        if (line.empty() || line[0] == '#') continue;

        size_t col = line.find(':');
        if (col == string::npos) {
            cerr << "Syntaxe invalide ligne " << nline << endl; return false;
        }

        int cycle = stoi(line.substr(0, col));
        string name = StringUtils::trim(line.substr(col + 1));
        if (cycle < prev) {
            cerr << "Cycles non croissants ligne " << nline << endl; return false;
        }

        finishUntil(cycle);

        const Process *p = nullptr;
        for (const auto &pr : config.getProcesses())
            if (pr.name == name) { p = &pr; break; }

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
