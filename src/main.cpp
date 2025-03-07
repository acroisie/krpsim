#include "Parser.hpp"
#include "Config.hpp"
#include "ProcessManager.hpp"
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <filename> <delay>" << std::endl;
        return 1;
    }
    
    Parser parser(argv[1]);
    Config config;
    if (parser.parse(config)) {
        cout << "Stocks after parsing: " << endl;
        for (const Stock &stock : config.getStocks()) {
            cout << "Stock: " << stock.name << " " << stock.quantity << endl;
        }
        int delay = std::stoi(argv[2]);
        ProcessManager processManager(config, delay);
        processManager.runSimulation();
    }
    return 0;
}