#include "Parser.hpp"
#include "Lexer.hpp"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <unordered_map>

using namespace std;

Parser::Parser(const string &filename) : filename_(filename), optimizeFound_(false) {}

bool Parser::parse(Config &config) {
    ifstream file(filename_);
    if (!file) {
        cerr << "Error: could not open file " << filename_ << endl;
        return false;
    }
    string line;
    int lineNumber = 0;

    try {
        while (getline(file, line)) {
            ++lineNumber;
            size_t pos = line.find_first_not_of(" \t");
            if (pos == string::npos || line[pos] == '#') {
                continue;
            }
            parseLine(line, config);
        }
        if (config.getProcesses().empty()) {
            throw runtime_error("No process defined in file.");
        }
    } catch (const exception &e) {
        cerr << "Error at line " << lineNumber << ": " << e.what() << endl;
        return false;
    }
    return true;
}

void Parser::parseLine(const string &line, Config &config) {
    Lexer lex(line);
    string firstToken = lex.nextIdentifier();

    if (firstToken == "optimize") {
        if (optimizeFound_) {
            throw runtime_error("Multiple optimize lines found.");
        }
        optimizeFound_ = true;
        lex.expect(':');
        lex.expect('(');
        auto optList = parseOptimizeList(lex);
        lex.expect(')');
        config.setOptimizeGoal(optList);
        return;
    }

    lex.expect(':');
    if (lex.peek() == '(') {
        Process proc;
        proc.name = firstToken;

        lex.expect('(');
        proc.inputs = parseResourceList(lex);
        lex.expect(')');
        lex.expect(':');

        lex.expect('(');
        proc.outputs = parseResourceList(lex);
        lex.expect(')');
        lex.expect(':');

        proc.nbCycle = lex.nextInteger();
        config.addProcess(proc);
    } else {
        Stock stock;
        stock.name = firstToken;
        stock.quantity = lex.nextInteger();
        config.addStock(stock);
    }
}

unordered_map<string, int> Parser::parseResourceList(Lexer &lex) {
    unordered_map<string, int> resources;
    while (true) {
        string name = lex.nextIdentifier();
        lex.expect(':');
        int quantity = lex.nextInteger();
        resources[name] = quantity;
        if (!lex.match(';')) {
            break;
        }
    }
    return resources;
}

vector<string> Parser::parseOptimizeList(Lexer &lex) {
    vector<string> optimizeList;
    while (true) {
        optimizeList.push_back(lex.nextIdentifier());
        if (!lex.match(';')) {
            break;
        }
    }
    return optimizeList;
}
