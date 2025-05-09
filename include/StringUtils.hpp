#pragma once
#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

namespace StringUtils {
    inline std::string trim(const std::string &str) {
        const auto start = str.find_first_not_of(" \t");
        if (start == std::string::npos) {
            return "";
        }
        const auto end = str.find_last_not_of(" \t");
        return str.substr(start, end - start + 1);
    }

    inline std::vector<std::string> split(const std::string &str,
                                          char delimiter) {
        std::vector<std::string> elements;
        std::istringstream stream(str);
        std::string element;
        while (std::getline(stream, element, delimiter)) {
            elements.push_back(element);
        }
        return elements;
    }
}