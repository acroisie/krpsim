#pragma once
#include <string>

class Lexer {
  public:
    explicit Lexer(const std::string &input);

    void skipWhitespace();
    char peek();
    void expect(char c);
    bool match(char c);
    std::string nextIdentifier();
    int nextInteger();

  private:
    const std::string &input;
    size_t position;
};
