#pragma once

#include "Parser.hpp"
#include "StringUtils.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>

using namespace std;

Parser::Parser(const string &filename) : _filename(filename) {}

bool Parser::parse()
{
	ifstream file{_filename};
	if (!file)
	{
		cerr << "Error: could not open file " << _filename << endl;
		return false;
	}

	string line;

	while (getline(file, line))
	{
		line = trim(line);
		if (line.empty() || line.front() == '#')
		{
			continue;
		}

		if (line.find("optimize:") != string::npos)
		{
			parseOptimization(line);
		}
		else if (line.find("stock:") != string::npos)
		{
			parseStock(line);
		}
		else if (line.find("process:") != string::npos)
		{
			parseProcess(line);
		}
		else
		{
			cerr << "Error: unknown token " << line << endl;
			return false;
		}
	}
	return true;
}

const vector<Stock> &Parser::getStocks() const
{
	return _stocks;
}

const vector<Process> &Parser::getProcesses() const
{
	return _processes;
}

void Parser::parseStock(const string& line) {
	auto elements = split(line, ':');
	if (elements.size() == 2) {
		cout << "Stock: " << elements[1] << endl; // TODO: remove this line
		_stocks.push_back(elements)
	}
}
