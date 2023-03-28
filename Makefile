CXX = g++
CXXFLAGS = -Wall -Werror -Wextra -pedantic -std=c++17 -pthread -O0
LDFLAGS = 

SRC = main.cpp
OBJ = $(SRC:.cc=.o)
EXEC = main

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $(OBJ) $(LBLIBS)

clean:
	rm -rf $(OBJ) $(EXEC)