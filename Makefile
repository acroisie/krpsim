NAME = krpsim
CXX = g++
CXXFLAGS = -Wall -Wextra -Werror -std=c++17 -Iinclude -g

SRC := $(shell find src -type f -name "*.cpp" ! -name "*tests*")

OBJS := $(SRC:.cpp=.o)

all: $(NAME)

$(NAME): $(OBJS)
	$(CXX) $(OBJS) -o $(NAME)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS)

fclean: clean
	rm -f $(NAME)

re: fclean all

.PHONY: all clean fclean re
