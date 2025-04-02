NAME = krpsim
CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++17 -Iinclude -O3

SRC_DIR = src
INCLUDE_DIR = include
OBJ_DIR = obj

SRCS = $(SRC_DIR)/main.cpp \
       $(SRC_DIR)/Config.cpp \
       $(SRC_DIR)/Lexer.cpp \
       $(SRC_DIR)/Parser.cpp \
       $(SRC_DIR)/Simulator.cpp \
       $(SRC_DIR)/GeneticAlgorithm.cpp \
       $(SRC_DIR)/ProcessManager.cpp

OBJS = $(SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

all: $(NAME)

$(NAME): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(NAME)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

clean:
	rm -rf $(OBJ_DIR)

fclean: clean
	rm -f $(NAME)

re: fclean all

.PHONY: all clean fclean re