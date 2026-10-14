# Minimal build file (no external dependencies).  `make test` runs all tests.
CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -g
DEPFLAGS := -MMD -MP
BUILD    := build

LIB_SRCS := tree.cpp single_tree.cpp naive_distance.cpp hdt.cpp count.cpp
LIB_OBJS := $(LIB_SRCS:%.cpp=$(BUILD)/%.o)

TEST_SRCS := $(wildcard tests/test_*.cpp)
TEST_BINS := $(TEST_SRCS:tests/%.cpp=$(BUILD)/%)

all: $(BUILD)/soda13 $(TEST_BINS)

$(BUILD)/%.o: %.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -c $< -o $@

$(BUILD)/soda13: main.cpp $(LIB_OBJS) | $(BUILD)
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -MF $@.d $< $(LIB_OBJS) -o $@

$(BUILD)/test_%: tests/test_%.cpp $(LIB_OBJS) | $(BUILD)
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -MF $@.d $< $(LIB_OBJS) -o $@

$(BUILD):
	mkdir -p $(BUILD)

test: $(TEST_BINS) $(BUILD)/soda13
	@for t in $(TEST_BINS); do echo "== $$t"; $$t || exit 1; done
	@echo "== tests/test_cli.sh"; sh tests/test_cli.sh $(BUILD)/soda13

clean:
	rm -rf $(BUILD)

-include $(wildcard $(BUILD)/*.d)

.PHONY: all test clean
