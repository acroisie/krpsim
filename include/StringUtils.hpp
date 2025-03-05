#pragma once
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>

inline std::string trim(const std::string &str) {
	const auto strBegin = str.find_first_not_of(" \t");
	if (strBegin == std::string::npos) {
		return "";
	}
	const auto strEnd = str.find_last_not_of(" \t");

	return str.substr(strBegin, strEnd - strBegin + 1);
}

inline std::vector<std::string> split(const std::string &str, char delimiter) {
	std::vector<std::string> elements;
	std::istringstream stream(str);
	std::string element;
	
	while (std::getline(stream, element, delimiter)) {
		elements.push_back(element);
	}

	return elements;
}