#include "Lexer.hpp"
#include <sstream>
#include <stdexcept>
#include <cctype>

Lexer::Lexer(const std::string& input) : input_(input), pos_(0) {}

void Lexer::skipWhitespace() {
    while (pos_ < input_.size() && std::isspace(input_[pos_]))
        ++pos_;
}

char Lexer::peek() {
    skipWhitespace();
    return (pos_ < input_.size()) ? input_[pos_] : '\0';
}

void Lexer::expect(char c) {
    skipWhitespace();
    if (peek() != c) {
        std::ostringstream oss;
        oss << "Expected '" << c << "' at position " << pos_;
        throw std::runtime_error(oss.str());
    }
    ++pos_;
}

bool Lexer::match(char c) {
    skipWhitespace();
    if (peek() == c) { ++pos_; return true; }
    return false;
}

std::string Lexer::nextIdentifier() {
    skipWhitespace();
    size_t start = pos_;
    while (pos_ < input_.size() && (std::isalnum(input_[pos_]) || input_[pos_] == '_' || input_[pos_] == '-'))
        ++pos_;
    if (start == pos_)
        throw std::runtime_error("Expecting identifier at position " + std::to_string(pos_));
    return input_.substr(start, pos_ - start);
}

int Lexer::nextInteger() {
    skipWhitespace();
    size_t start = pos_;
    if (pos_ >= input_.size() || !std::isdigit(input_[pos_]))
        throw std::runtime_error("Expecting integer at position " + std::to_string(pos_));
    while (pos_ < input_.size() && std::isdigit(input_[pos_]))
        ++pos_;
    if (start == pos_)
        throw std::runtime_error("Expecting integer at position " + std::to_string(pos_));
    try {
        return std::stoi(input_.substr(start, pos_ - start));
    } catch (...) {
        throw std::runtime_error("Invalid format for integer at position " + std::to_string(pos_));
    }
}