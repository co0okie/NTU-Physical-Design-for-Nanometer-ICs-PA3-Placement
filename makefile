.PHONY: all clean

CXX=g++
CXXFLAGS=-std=c++17 -static -O2 -Wall -D_GLIBCXX_ISE_CXX11_ABI=1
LDFLAGS=-Llib -lDetailPlace -lGlobalPlace -lLegalizer -lPlacement -lParser -lPlaceCommon
SRCS=src/ObjectiveFunction.cpp src/Optimizer.cpp src/GlobalPlacer.cpp src/main.cpp
OBJS=$(SRCS:.cpp=.o)
DEPS=$(OBJS:.o=.d)
INCLUDES=$(wildcard src/*.h)
EXE=bin/place

all: $(EXE)

$(EXE): $(OBJS)
	mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

-include $(DEPS)

clean:
	rm -rf src/*.o src/*.d bin