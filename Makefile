# Builds sesmesh into bin/ with a C++17 compiler.  It runs in parallel on std::threads, so it needs
# nothing beyond the C++ standard library.
#
#   make                  bin/sesmesh
#   make CXX=clang++      another compiler
#   make test             build and run the unit tests

CXX      = g++
CXXFLAGS = -O2 -std=c++17 -Wall -Wextra -pthread

SRC = src/ses.cpp src/mesh.cpp src/clean.cpp src/io.cpp src/parallel.cpp
HDR = src/vec3.h src/ses.h src/mesh.h src/clean.h src/io.h src/parallel.h

all: bin/sesmesh

bin:
	mkdir -p bin

bin/sesmesh: src/sesmesh.cpp $(SRC) $(HDR) | bin
	$(CXX) $(CXXFLAGS) -o $@ src/sesmesh.cpp $(SRC)

bin/test_sesmesh: test/test_sesmesh.cpp $(SRC) $(HDR) | bin
	$(CXX) $(CXXFLAGS) -Isrc -o $@ test/test_sesmesh.cpp $(SRC)

test: bin/test_sesmesh
	bin/test_sesmesh

clean:
	rm -rf bin

.PHONY: all test clean
