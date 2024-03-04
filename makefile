BIN := bin
SRC := synthetiq
INC := include/eigen-3.3.9/
OBJ := build

app     := $(BIN)/main $(BIN)/comparison_generator $(BIN)/main_resynth
sources := $(wildcard $(SRC)/*.h)
objects := $(subst $(SRC),$(OBJ),$(sources:.h=.o))
sources_all := $(wildcard $(SRC)/*.cpp)
deps    := $(sources_all:.cpp=.d)

CXX      := g++
CPPFLAGS := -I $(INC) -MMD -MP -march=native -Ofast -ffast-math -DEIGEN_NO_DEBUG -funroll-loops -fprefetch-loop-arrays -mtune=native -flto=6 -frename-registers
CXXFLAGS := -std=c++17 -fopenmp

.PHONY: all directories clean

directories:
	mkdir -p $(BIN) #
	mkdir -p $(OBJ) #

all: directories $(app)

$(BIN)/% : $(objects) $(OBJ)/%.o
	$(CXX) $(CXXFLAGS) $^ -o $@

$(OBJ)/%.o: $(SRC)/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $^ -o $@

clean:
	@$(RM) $(app) $(OBJ)/*.d $(OBJ)/*.o

-include $(deps)