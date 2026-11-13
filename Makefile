# Builds sesmesh into bin/ with a C++17 compiler.
#
#   make                  with OpenMP
#   make OPENMP=          without OpenMP (one thread; macOS clang has none by default)
#   make CXX=clang++      another compiler
#   make test             build and run the unit tests

CXX      = g++
CXXFLAGS = -O2 -std=c++17 -Wall -Wextra
OPENMP   = -fopenmp
OMPFLAGS = $(if $(OPENMP),$(OPENMP),-Wno-unknown-pragmas)

SRC = src/ses.cpp src/mesh.cpp src/clean.cpp src/io.cpp
HDR = src/vec3.h src/ses.h src/mesh.h src/clean.h src/io.h

all: bin/sesmesh

bin:
	mkdir -p bin

bin/sesmesh: src/sesmesh.cpp $(SRC) $(HDR) | bin
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) -o $@ src/sesmesh.cpp $(SRC)

bin/test_sesmesh: test/test_sesmesh.cpp $(SRC) $(HDR) | bin
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) -Isrc -o $@ test/test_sesmesh.cpp $(SRC)

test: bin/test_sesmesh
	bin/test_sesmesh

clean:
	rm -rf bin

.PHONY: all test clean
