CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2 -I./include
LDFLAGS = 

# Directories
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin

# Source files
SRC_COMMON = $(wildcard $(SRC_DIR)/*.cpp)
SRC_COMMON := $(filter-out $(SRC_DIR)/krpsim.cpp $(SRC_DIR)/krpsim_verif.cpp, $(SRC_COMMON))

# Object files
OBJ_COMMON = $(SRC_COMMON:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
OBJ_KRPSIM = $(OBJ_COMMON) $(OBJ_DIR)/krpsim.o
OBJ_VERIF = $(OBJ_COMMON) $(OBJ_DIR)/krpsim_verif.o

# Binaries
BIN_KRPSIM = $(BIN_DIR)/krpsim
BIN_VERIF = $(BIN_DIR)/krpsim_verif

# Targets
all: $(BIN_KRPSIM) $(BIN_VERIF)

# Create bin directory
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# Create obj directory
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Compile krpsim
$(BIN_KRPSIM): $(OBJ_KRPSIM) | $(BIN_DIR)
	$(CXX) $(LDFLAGS) $^ -o $@

# Compile krpsim_verif
$(BIN_VERIF): $(OBJ_VERIF) | $(BIN_DIR)
	$(CXX) $(LDFLAGS) $^ -o $@

# Compile common source files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean objects
clean:
	rm -rf $(OBJ_DIR)

# Clean all
fclean: clean
	rm -rf $(BIN_DIR)

# Rebuild all
re: fclean all

.PHONY: all clean fclean re