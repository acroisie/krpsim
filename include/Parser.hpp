#pragma once
#include <string>
#include <vector>
#include "Config.hpp"

class Parser {
public:
    explicit Parser(const std::string &filename);
    bool parse(class Config &config);

private:
    std::string filename_;
    bool optimizeFound_;

    void parseLine(const std::string &line, class Config &config);
    std::unordered_map<std::string, int> parseResourceList(class Lexer &lex);
    std::vector<std::string> parseOptimizeList(class Lexer &lex);
};
