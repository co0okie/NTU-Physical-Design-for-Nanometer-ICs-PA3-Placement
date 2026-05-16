.PHONY: all clean run

CXX=g++
CXXFLAGS=-std=c++17 -static -O3 -Wall -MMD -MP -D_GLIBCXX_ISE_CXX11_ABI=1
# CXXFLAGS=-std=c++17 -static -O0 -g -Wall -MMD -MP -D_GLIBCXX_ISE_CXX11_ABI=1
LDFLAGS=-Llib -lDetailPlace -lGlobalPlace -lLegalizer -lPlacement -lParser -lPlaceCommon
SRCS=src/ObjectiveFunction.cpp src/Optimizer.cpp src/GlobalPlacer.cpp src/main.cpp
OBJS=$(SRCS:.cpp=.o)
DEPS=$(OBJS:.o=.d)
INCLUDES=$(wildcard src/*.h)
EXE=bin/place

# BENCHMARK=benchmark/ibm01/ibm01-cu85.aux
BENCHMARK=benchmark/ibm05/ibm05.aux

TGZ=b11107051_pa3.tgz

all: $(EXE)

run: $(EXE)
	./$< -aux $(BENCHMARK)

tgz: $(TGZ)

%.tgz: $(EXE) $(SRCS) $(INCLUDES) makefile readme.txt report/report.pdf
	temp_dir=$$(mktemp -d); \
	mkdir -p "$$temp_dir"/$*/src "$$temp_dir"/$*/bin "$$temp_dir"/$*/lib; \
	cp src/*.h src/*.cpp "$$temp_dir"/$*/src/; \
	cp lib/* "$$temp_dir"/$*/lib/; \
	cp $(EXE) "$$temp_dir"/$*/bin/; \
	cp makefile readme.txt report/report.pdf "$$temp_dir"/$*; \
	tar zcvf $@ -C "$$temp_dir" $*; \
	rm -rf "$$temp_dir"

$(EXE): $(OBJS)
	mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

-include $(DEPS)

clean:
	rm -rf src/*.o src/*.d bin/ *.dp.pl *.gp.pl *.lg.pl init.plt $(TGZ) \
	report/report.aux report/report.fdb_latexmk report/report.fls report/report.log report/report.out report/report.synctex.gz report/report.xdv