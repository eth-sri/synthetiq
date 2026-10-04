# Keep the original C++ entry points; install Rust alongside them.
.DEFAULT_GOAL := all
CARGO ?= cargo
CXX ?= g++
BIN := bin
OBJ := build/cpp
SRC := synthetiq
CPPFLAGS := -I include/eigen-3.3.9/ -MMD -MP -march=native -Ofast -ffast-math -DEIGEN_NO_DEBUG -funroll-loops -fprefetch-loop-arrays -mtune=native -flto=6 -frename-registers
CXXFLAGS := -std=c++17 -fopenmp
headers := $(wildcard $(SRC)/*.h)
objects := $(patsubst $(SRC)/%.h,$(OBJ)/%.o,$(headers))
cpp_apps := $(BIN)/cpp/main $(BIN)/cpp/main_resynth $(BIN)/cpp/comparison_generator

.PHONY: all rust cpp test clean
all: cpp
rust:
	$(CARGO) build --release --bin synthetiq --bin main_resynth --bin comparison_generator
	mkdir -p $(BIN)
	cp target/release/synthetiq $(BIN)/rust
	cp target/release/main_resynth $(BIN)/rust_resynth
	cp target/release/comparison_generator $(BIN)/rust_comparison_generator

cpp: $(cpp_apps)
	cp $(BIN)/cpp/main $(BIN)/main
	cp $(BIN)/cpp/main_resynth $(BIN)/main_resynth
	cp $(BIN)/cpp/comparison_generator $(BIN)/comparison_generator
$(BIN)/cpp/%: $(objects) $(OBJ)/%.o
	mkdir -p $(BIN)/cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

$(OBJ)/%.o: $(SRC)/%.cpp
	mkdir -p $(OBJ)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

-include $(wildcard $(OBJ)/*.d)
test:
	$(CARGO) test --all-targets
clean:
	$(RM) $(BIN)/main $(BIN)/main_resynth $(BIN)/comparison_generator $(BIN)/rust $(BIN)/rust_resynth $(BIN)/rust_comparison_generator $(cpp_apps)
	$(RM) -r target build/cpp build/reference
