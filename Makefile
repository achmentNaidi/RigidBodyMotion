CXX       := g++ -ggdb
CXX_FLAGS := -std=c++17

BIN     := bin
SRC     := src
INCLUDE := include

LIBRARIES   :=
EXECUTABLE  := main

RED = \033[1;31m
GREEN = \033[1;32m
BLUE = \033[1;34m
YELLOW = \033[1;33m
NC = \033[1;0m

all: $(BIN)/$(EXECUTABLE)

run: clean all
	clear
	@echo "$(YELLOW)Running...$(NC)"
	./$(BIN)/$(EXECUTABLE)

$(BIN)/$(EXECUTABLE): $(SRC)/*.cpp
	@echo "$(GREEN)Build...$(NC)"
	$(CXX) $(CXX_FLAGS) -I$(INCLUDE) $^ -o $@ $(LIBRARIES)
	@echo "$(RED)Finished!$(NC)"

clean:
	-rm $(BIN)/*