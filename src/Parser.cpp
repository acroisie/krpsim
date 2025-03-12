#include "Parser.hpp"
#include "Lexer.hpp"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <unordered_map>

Parser::Parser(const std::string &filename)
    : filename_(filename), optimizeFound_(false) {}

bool Parser::parse(Config &config) {
    std::ifstream file(filename_);
    if (!file) {
        std::cerr << "Error: could not open file " << filename_ << std::endl;
        return false;
    }
    std::string line;
    int lineNumber = 0;

    try {
        while (std::getline(file, line)) {
            ++lineNumber;
            // Ignore lignes vides ou commentaires
            size_t pos = line.find_first_not_of(" \t");
            if (pos == std::string::npos || line[pos] == '#') {
                continue;
            }
            parseLine(line, config);
        }
        if (config.getProcesses().empty()) {
            throw std::runtime_error("No process defined in file.");
        }
    } catch (const std::exception &e) {
        std::cerr << "Error at line " << lineNumber << ": " << e.what() << std::endl;
        return false;
    }
    return true;
}

void Parser::parseLine(const std::string &line, Config &config) {
    Lexer lex(line);
    std::string firstToken = lex.nextIdentifier();

    if (firstToken == "optimize") {
        if (optimizeFound_) {
            throw std::runtime_error("Multiple optimize lines found.");
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
        // C’est un process
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
        // C’est un stock
        Stock stock;
        stock.name = firstToken;
        stock.quantity = lex.nextInteger();
        config.addStock(stock);
    }
}

std::unordered_map<std::string, int> Parser::parseResourceList(Lexer &lex) {
    std::unordered_map<std::string, int> resources;
    while (true) {
        std::string name = lex.nextIdentifier();
        lex.expect(':');
        int quantity = lex.nextInteger();
        resources[name] = quantity;
        if (!lex.match(';')) {
            break;
        }
    }
    return resources;
}

std::vector<std::string> Parser::parseOptimizeList(Lexer &lex) {
    std::vector<std::string> optimizeList;
    while (true) {
        optimizeList.push_back(lex.nextIdentifier());
        if (!lex.match(';')) {
            break;
        }
    }
    return optimizeList;
}
