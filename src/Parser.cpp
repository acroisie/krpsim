#include "Parser.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <cctype>

using namespace std;

class Lexer {
public:
    explicit Lexer(const string& input) : input_(input), pos_(0) {}

    void skipWhitespace() {
        while (pos_ < input_.size() && isspace(input_[pos_]))
            ++pos_;
    }

    char peek() {
        skipWhitespace();
        return (pos_ < input_.size()) ? input_[pos_] : '\0';
    }

    void expect(char c) {
        skipWhitespace();
        if (peek() != c) {
            ostringstream oss;
            oss << "Expected '" << c << "' at position " << pos_;
            throw runtime_error(oss.str());
        }
        ++pos_;
    }

    bool match(char c) {
        skipWhitespace();
        if (peek() == c) { ++pos_; return true; }
        return false;
    }

    string nextIdentifier() {
        skipWhitespace();
        size_t start = pos_;
        while (pos_ < input_.size() && (isalnum(input_[pos_]) || input_[pos_] == '_' || input_[pos_] == '-'))
            ++pos_;
        if (start == pos_)
            throw runtime_error("Expecting identifier at position " + to_string(pos_));
        return input_.substr(start, pos_ - start);
    }

    int nextInteger() {
        skipWhitespace();
        size_t start = pos_;
        if (!isdigit(input_[pos_]))
            throw runtime_error("Expecting integer at position " + to_string(pos_));
        while (pos_ < input_.size() && isdigit(input_[pos_]))
            ++pos_;
        if (start == pos_)
            throw runtime_error("Expecting integer at position " + to_string(pos_));
        try {
            return stoi(input_.substr(start, pos_ - start));
        } catch (...) {
            throw runtime_error("Invalid format for integer at position " + to_string(pos_));
        }
    }

private:
    const string& input_;
    size_t pos_;
};

static unordered_map<string, int> parseResourceList(Lexer& lex) {
    unordered_map<string, int> resources;
    while (true) {
        string name = lex.nextIdentifier();
        lex.expect(':');
        int quantity = lex.nextInteger();
        resources[name] = quantity;
        if (!lex.match(';'))
            break;
    }
    return resources;
}

static vector<string> parseOptimizeList(Lexer& lex) {
    vector<string> optimizeList;
    while (true) {
        optimizeList.push_back(lex.nextIdentifier());
        if (!lex.match(';'))
            break;
    }
    return optimizeList;
}

Parser::Parser(const string& filename) : filename_(filename), optimizeFound_(false) {}

bool Parser::parse() {
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
            string trimmed = line;
            size_t pos = trimmed.find_first_not_of(" \t");
            if (pos == string::npos || trimmed[pos] == '#')
                continue;
            parseLine(trimmed, lineNumber);
        }
        if (processes_.empty())
            throw runtime_error("No process defined in file.");
    } catch (const exception& e) {
        cerr << "Error at line " << lineNumber << ": " << e.what() << endl;
        return false;
    }
    return true;
}

const vector<Stock>& Parser::getStocks() const {
    return stocks_;
}

const vector<Process>& Parser::getProcesses() const {
    return processes_;
}

void Parser::parseLine(const string& line, int lineNumber) {
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
        cout << "Optimize: " << lineNumber << endl;
        for (const auto& token : optList)
            cout << token << " ";
        cout << endl;
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
        processes_.push_back(proc);
    } else {
        Stock stock;
        stock.name = firstToken;
        stock.quantity = lex.nextInteger();
        stocks_.push_back(stock);
    }
}
