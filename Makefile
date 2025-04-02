NAME = krpsim
CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++17 -Iinclude -O3 -DNDEBUG

SRC := $(shell find src -type f -name "*.cpp" ! -name "*tests*")

OBJS := $(SRC:.cpp=.o)

all: $(NAME)

$(NAME): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(NAME)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS)

fclean: clean
	rm -f $(NAME)

re: fclean all

.PHONY: all clean fclean re
