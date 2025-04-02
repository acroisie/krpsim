#include "Lexer.hpp"
#include <cctype>
#include <stdexcept>

using namespace std;

Lexer::Lexer(const string &input) : input(input), position(0) {}

void Lexer::skipWhitespace() {
    while (position < input.size() &&
           isspace(static_cast<unsigned char>(input[position]))) {
        ++position;
    }
}

char Lexer::peek() {
    skipWhitespace();
    if (position < input.size()) {
        return input[position];
    }
    return '\0';
}

void Lexer::expect(char c) {
    skipWhitespace();
    if (peek() != c) {
        throw runtime_error("Expected '" + string(1, c) + "' at position " +
                            to_string(position));
    }
    ++position;
}

bool Lexer::match(char c) {
    skipWhitespace();
    if (peek() == c) {
        ++position;
        return true;
    }
    return false;
}

string Lexer::nextIdentifier() {
    skipWhitespace();
    size_t start = position;
    while (position < input.size() &&
           (isalnum(static_cast<unsigned char>(input[position])) ||
            input[position] == '_' || input[position] == '-')) {
        ++position;
    }
    if (start == position) {
        throw runtime_error("Expecting identifier at position " +
                            to_string(position));
    }
    return input.substr(start, position - start);
}

int Lexer::nextInteger() {
    skipWhitespace();
    size_t start = position;
    if (position >= input.size() ||
        !isdigit(static_cast<unsigned char>(input[position]))) {
        throw runtime_error("Expecting integer at position " +
                            to_string(position));
    }
    while (position < input.size() &&
           isdigit(static_cast<unsigned char>(input[position]))) {
        ++position;
    }
    return stoi(input.substr(start, position - start));
}