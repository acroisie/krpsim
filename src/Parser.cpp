#include "../include/parser.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>

using namespace std;

Parser::Parser(string &filename)
{
	ifstream file(filename);
	if (!file.is_open())
	{
		throw runtime_error("Cannot open file:" + filename);
	}

	string line;
	while (getline(file, line))
	{
		if (line.empty() || line[0] == '#') continue;
		
	}
}