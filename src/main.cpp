#include "Parser.hpp"
#include "Config.hpp"
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <filename>" << " " << "<delay>"<< endl;
        return 1;
    }
    
    Parser parser(argv[1]);
    Config config;
    if (parser.parse(config)) {
        cout << endl << "Stocks :" << endl;
        for (const auto& stock : config.getStocks()) {
            cout << "  " << stock.name << " => " << stock.quantity << endl;
        }
        cout << endl << "Process :" << endl;
        for (const auto& proc : config.getProcesses()) {
            cout << "  " << proc.name << " (cycle: " << proc.nbCycle << ")" << endl;
            cout << "    Inputs:" << endl;
            for (const auto& [res, qty] : proc.inputs)
                cout << "      " << res << " : " << qty << endl;
            cout << "    Outputs:" << endl;
            for (const auto& [res, qty] : proc.outputs)
                cout << "      " << res << " : " << qty << endl;
        }
         cout << endl << "Optimize Goal :" << endl;
        for (const auto& optimize : config.getOptimizeGoal()) {
            cout << "  " << optimize << endl;
        }
    }
    return 0;
}