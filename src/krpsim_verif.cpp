#include "Parser.hpp"
#include "TraceVerifier.hpp"
#include <iostream>

int main(int argc, char *argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <config_file> <trace_file>\n";
        return 1;
    }

    Parser parser(argv[1]);
    Config config;
    if (!parser.parse(config)) return 1;

    TraceVerifier verifier(config);          // limite par défaut 1 000 000 cycles
    bool ok = verifier.verifyFile(argv[2]);
    return ok ? 0 : 1;
}
