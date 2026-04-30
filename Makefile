# Compiler and flags
CXX      := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -O2

# Directories
SRC_DIR := src
BIN_DIR := bin

# Target
TARGET := $(BIN_DIR)/sketch_test

# Sources
SRCS := $(SRC_DIR)/main.cpp $(SRC_DIR)/models.cpp

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRCS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $(SRCS)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

clean:
	rm -f $(TARGET)
	rm -f $(BIN_DIR)/*.csv
