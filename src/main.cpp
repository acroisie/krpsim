#include "Config.hpp"
#include "Parser.hpp"
#include "ProcessManager.hpp"
#include <cstdlib>
#include <iostream>
using namespace std;

int main(int argc, char *argv[]) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <filename> <delay>" << endl;
        return 1;
    }

    Parser parser(argv[1]);
    Config config;

    if (!parser.parse(config)) {
        return 1;
    }

    int delay = atoi(argv[2]);
    if (delay <= 0) {
        cerr << "Error: Delay must be a positive integer." << endl;
        return 1;
    }

    ProcessManager manager(config, delay);
    manager.run();

    return 0;
}