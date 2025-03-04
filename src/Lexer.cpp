#include "Lexer.hpp"
#include <sstream>
#include <stdexcept>
#include <cctype>

using namespace std;

Lexer::Lexer(const string& input) : input_(input), pos_(0) {}

void Lexer::skipWhitespace() {
    while (pos_ < input_.size() && isspace(input_[pos_]))
        ++pos_;
}

char Lexer::peek() {
    skipWhitespace();
    return (pos_ < input_.size()) ? input_[pos_] : '\0';
}

void Lexer::expect(char c) {
    skipWhitespace();
    if (peek() != c) {
        ostringstream oss;
        oss << "Expected '" << c << "' at position " << pos_;
        throw runtime_error(oss.str());
    }
    ++pos_;
}

bool Lexer::match(char c) {
    skipWhitespace();
    if (peek() == c) { ++pos_; return true; }
    return false;
}

string Lexer::nextIdentifier() {
    skipWhitespace();
    size_t start = pos_;
    while (pos_ < input_.size() && (isalnum(input_[pos_]) || input_[pos_] == '_' || input_[pos_] == '-'))
        ++pos_;
    if (start == pos_)
        throw runtime_error("Expecting identifier at position " + to_string(pos_));
    return input_.substr(start, pos_ - start);
}

int Lexer::nextInteger() {
    skipWhitespace();
    size_t start = pos_;
    if (pos_ >= input_.size() || !isdigit(input_[pos_]))
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