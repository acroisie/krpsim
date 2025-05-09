CXX      = g++
CXXFLAGS = -Wall -Wextra -std=c++17 -Iinclude -O3

SRC_DIR  = src
OBJ_DIR  = obj

COMMON_SRCS = \
	$(SRC_DIR)/Config.cpp \
	$(SRC_DIR)/Lexer.cpp  \
	$(SRC_DIR)/Parser.cpp \
	$(SRC_DIR)/Simulator.cpp

KRPSIM_SRCS = \
	$(SRC_DIR)/main.cpp \
	$(SRC_DIR)/GeneticAlgorithm.cpp \
	$(SRC_DIR)/ProcessManager.cpp

KRPSIM_OBJS = $(COMMON_SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o) \
	          $(KRPSIM_SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

VERIF_SRCS  = \
	$(SRC_DIR)/TraceVerifier.cpp \
	$(SRC_DIR)/krpsim_verif.cpp

VERIF_OBJS  = $(COMMON_SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o) \
	          $(VERIF_SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

all: krpsim krpsim_verif

krpsim: $(KRPSIM_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

krpsim_verif: $(VERIF_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

clean:
	rm -rf $(OBJ_DIR)

fclean: clean
	rm -f krpsim krpsim_verif tracefile.txt

re: fclean all

.PHONY: all clean fclean re
