#include "Parser.hpp"
#include <iostream>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage : " << argv[0] << " <file>\n";
        return 1;
    }
    Parser parser{argv[1]};
    if (parser.parse()) {
        std::cout << "Stocks lus :\n";
        for (const auto& stock : parser.getStocks())
            std::cout << "  " << stock.name << " => " << stock.quantity << "\n";
        std::cout << "Processus lus :\n";
        for (const auto& proc : parser.getProcesses()) {
            std::cout << "  " << proc.name << " (cycle: " << proc.cycleAmount << ")\n";
            std::cout << "    Inputs:\n";
            for (const auto& [res, qty] : proc.inputs)
                std::cout << "      " << res << " : " << qty << "\n";
            std::cout << "    Outputs:\n";
            for (const auto& [res, qty] : proc.outputs)
                std::cout << "      " << res << " : " << qty << "\n";
        }
    }
    return 0;
}
