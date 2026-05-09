.PHONY: all clean run

CXX=g++
CXXFLAGS=-std=c++17 -static -O2 -Wall -MMD -MP -D_GLIBCXX_ISE_CXX11_ABI=1
LDFLAGS=-Llib -lDetailPlace -lGlobalPlace -lLegalizer -lPlacement -lParser -lPlaceCommon
SRCS=src/ObjectiveFunction.cpp src/Optimizer.cpp src/GlobalPlacer.cpp src/main.cpp
OBJS=$(SRCS:.cpp=.o)
DEPS=$(OBJS:.o=.d)
INCLUDES=$(wildcard src/*.h)
EXE=bin/place

# BENCHMARK=benchmark/ibm01/ibm01-cu85.aux
BENCHMARK=benchmark/ibm05/ibm05.aux

all: $(EXE)

run: $(EXE)
	./$< -aux benchmark/ibm01/ibm01-cu85.aux

$(EXE): $(OBJS)
	mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

-include $(DEPS)

clean:
	rm -rf src/*.o src/*.d bin