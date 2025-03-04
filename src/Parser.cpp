#include "Parser.hpp"
#include "Stock.hpp"
#include "Process.hpp"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <sstream>

class LineParser
{
public:
	explicit LineParser(const std::string &input) : _input(input), _pos(0) {}

	void skipWhithespace()
	{
		while (_pos < _input.size() && std::isspace(_input[_pos]))
		{
			++_pos;
		}
	}

	char peek()
	{
		skipWhithespace();
		return _pos < _input.size() ? _input[_pos] : '\0';
	}

	void expect(char c)
	{
		skipWhithespace();
		if (peek() != c)
		{
			throw std::runtime_error("Unexpected character");
		}
		++_pos;
	}

	bool match(char c)
	{
		skipWhithespace();
		if (peek() == c)
		{
			++_pos;
			return true;
		}
		return false;
	}

	std::string parseIdentifier()
	{
		skipWhithespace();
		size_t start = _pos;
		while (_pos < _input.size() && (std::isalnum(_input[_pos]) || _input[_pos] == '_' || _input[_pos] == '-'))
		{
			++_pos;
		}
		if (start == _pos)
		{
			throw std::runtime_error("Expected identifier");
		}
		return _input.substr(start, _pos - start);
	}

	int parseInteger()
	{
		skipWhithespace();
		size_t start = _pos;
		if (_pos < _input.size() && (_input[_pos] == '+' || _input[_pos] == '-'))
		{
			++_pos;
		}
		while (_pos < _input.size() && std::isdigit(_input[_pos]))
		{
			++_pos;
		}
		if (start == _pos)
		{
			throw std::runtime_error("Expected integer");
		}
		try
		{
			return std::stoi(_input.substr(start, _pos - start));
		}
		catch (const std::exception &e)
		{
			std::cerr << e.what() << std::endl;
		}

		return std::stoi(_input.substr(start, _pos - start));
	}

	bool eof() const
	{
		return _pos >= _input.size();
	}

private:
	const std::string &_input;
	size_t _pos;
};

std::unordered_map<std::string, int> parseResourceList(LineParser &parser)
{
	std::unordered_map<std::string, int> resources;
	while (parser.eof() == false)
	{
		parser.skipWhithespace();
		if (parser.peek() == ')')
		{
			break;
		}
		std::string resource = parser.parseIdentifier();
		parser.skipWhithespace();
		parser.expect(':');
		parser.skipWhithespace();
		int amount = parser.parseInteger();
		resources[resource] = amount;
		parser.skipWhithespace();

		if (parser.match(';'))
		{
			break;
		}
	}
	return resources;
}

std::vector<std::string> parseOptimizeLine(LineParser &parser)
{
	std::vector<std::string> optimize;
	while (parser.eof() == false)
	{
		parser.skipWhithespace();
		if (parser.peek() == ')')
		{
			break;
		}
		optimize.push_back(parser.parseIdentifier());
		parser.skipWhithespace();
		if (parser.match(';'))
		{
			break;
		}
	}
	return optimize;
}

Parser::Parser(const std::string &filename) : _filename(filename) {}

bool Parser::parse()
{
	std::ifstream file(_filename);
	if (!file)
	{
		std::cerr << "Error: could not open file " << _filename << std::endl;
		return false;
	}
	std::string line;
	int lineNumber = 0;

	while (std::getline(file, line))
	{
		++lineNumber;
		std::istringstream stream(line);
		std::string trimmed;
		std::getline(stream, trimmed);
		std::string temp = trimmed;
		size_t pos = temp.find_first_not_of(" \t");
		if (pos == std::string::npos || temp[pos] == '#')
		{
			continue;
		}
		try
		{
			parseLine(trimmed, lineNumber);
		}
		catch (const std::exception &e)
		{
			std::cerr << "Error in line " << lineNumber << ": " << e.what() << std::endl;
			return false;
		}
	}
	return true;
}

const std::vector<Stock> &Parser::getStocks() const
{
	return _stocks;
}

const std::vector<Process> &Parser::getProcesses() const
{
	return _processes;
}

void Parser::parseLine(const std::string &line, int lineNumber)
{
	LineParser parser(line);
	parser.skipWhithespace();

	std::string firstToken = parser.parseIdentifier();

	if (firstToken == "optimize")
	{
		parser.skipWhithespace();
		parser.expect(':');
		parser.skipWhithespace();
		parser.expect('(');

		auto optimize = parseOptimizeLine(parser);

		parser.expect(')');
	}
	else
	{
		parser.skipWhithespace();
		parser.expect(':');
		parser.skipWhithespace();

		if (parser.peek() == '(')
		{
			Process process;
			process.name = firstToken;
			parser.expect('(');
			process.inputs = parseResourceList(parser);
			parser.expect(')');
			parser.skipWhithespace();
			parser.expect(':');
			parser.skipWhithespace();
			parser.expect('(');
			process.outputs = parseResourceList(parser);
			parser.expect(')');
			parser.skipWhithespace();
			parser.expect(':');
			parser.skipWhithespace();
			process.cycleAmount = parser.parseInteger();
			_processes.push_back(process);
		}
		else
		{
			Stock stock;
			stock.name = firstToken;
			parser.expect('(');
			stock.quantity = parser.parseInteger();
			parser.expect(')');
			_stocks.push_back(stock);
		}
	}