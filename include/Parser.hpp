#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include "Config.hpp"

class Parser {
public:
    explicit Parser(const std::string &filename);
    bool parse(Config &config);

private:
    std::string filename_;
    bool optimizeFound_;

    void parseLine(const std::string &line, Config &config);
    std::unordered_map<std::string, int> parseResourceList(class Lexer &lex);
    std::vector<std::string> parseOptimizeList(class Lexer &lex);
};
