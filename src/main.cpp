#include "Parser.hpp"
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <filename>" << endl;
        return 1;
    }
    
    Parser parser(argv[1]);
    if (parser.parse()) {
        cout << endl << "Stocks :" << endl;
        for (const auto& stock : parser.getStocks()) {
            cout << "  " << stock.name << " => " << stock.quantity << endl;
        }
        cout << endl << "Process :" << endl;
        for (const auto& proc : parser.getProcesses()) {
            cout << "  " << proc.name << " (cycle: " << proc.nbCycle << ")" << endl;
            cout << "    Inputs:" << endl;
            for (const auto& [res, qty] : proc.inputs)
                cout << "      " << res << " : " << qty << endl;
            cout << "    Outputs:" << endl;
            for (const auto& [res, qty] : proc.outputs)
                cout << "      " << res << " : " << qty << endl;
        }
    }
    return 0;
}
