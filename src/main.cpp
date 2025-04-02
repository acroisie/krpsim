#include "Config.hpp"
#include "Parser.hpp"
#include "ProcessManager.hpp"
#include <cstdlib>
#include <iostream>

int main(int argc, char *argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <filename> <delay>" << std::endl;
        return 1;
    }

    Parser parser(argv[1]);
    Config config;

    if (!parser.parse(config)) {
        return 1;
    }

    int delay = std::atoi(argv[2]);
    if (delay <= 0) {
        std::cerr << "Error: Delay must be a positive integer." << std::endl;
        return 1;
    }

    ProcessManager manager(config, delay);
    manager.run();

    return 0;
}