#include "Lexer.hpp"
#include <cctype>
#include <sstream>
#include <stdexcept>

Lexer::Lexer(const std::string &input) : input_(input), pos_(0) {}

void Lexer::skipWhitespace() {
    while (pos_ < input_.size() && isspace(static_cast<unsigned char>(input_[pos_]))) {
        ++pos_;
    }
}

char Lexer::peek() {
    skipWhitespace();
    if (pos_ < input_.size()) {
        return input_[pos_];
    }
    return '\0';
}

void Lexer::expect(char c) {
    skipWhitespace();
    if (peek() != c) {
        throw std::runtime_error("Expected '" + std::string(1, c) + "' at position " + std::to_string(pos_));
    }
    ++pos_;
}

bool Lexer::match(char c) {
    skipWhitespace();
    if (peek() == c) {
        ++pos_;
        return true;
    }
    return false;
}

std::string Lexer::nextIdentifier() {
    skipWhitespace();
    size_t start = pos_;
    while (pos_ < input_.size() &&
           (std::isalnum(static_cast<unsigned char>(input_[pos_])) || input_[pos_] == '_' || input_[pos_] == '-')) {
        ++pos_;
    }
    if (start == pos_) {
        throw std::runtime_error("Expecting identifier at position " + std::to_string(pos_));
    }
    return input_.substr(start, pos_ - start);
}

int Lexer::nextInteger() {
    skipWhitespace();
    size_t start = pos_;
    if (pos_ >= input_.size() || !std::isdigit(static_cast<unsigned char>(input_[pos_]))) {
        throw std::runtime_error("Expecting integer at position " + std::to_string(pos_));
    }
    while (pos_ < input_.size() && std::isdigit(static_cast<unsigned char>(input_[pos_]))) {
        ++pos_;
    }
    try {
        return std::stoi(input_.substr(start, pos_ - start));
    } catch (...) {
        throw std::runtime_error("Invalid integer format at position " + std::to_string(pos_));
    }
}
