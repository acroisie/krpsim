NAME = krpsim
CXX = g++
CXXFLAGS = -Wall -Wextra -Werror -std=c++17 -Iinclude
SRC = src/main.cpp src/Parser.cpp
OBJ = $(SRC:.cpp=.o)

all: $(NAME)

$(NAME): $(OBJ)
	$(CXX) $(OBJ) -o $(NAME)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ)

fclean: clean
	rm -f $(NAME)

re: fclean all
