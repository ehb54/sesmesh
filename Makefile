# Builds sesmesh and rcoal into bin/ with a C++17 compiler.  sesmesh runs in parallel on
# std::threads, so it needs nothing beyond the C++ standard library.
#
#   make                  bin/sesmesh, bin/rcoal and bin/rcoal3
#   make CXX=clang++      another compiler
#   make test             build and run the unit tests and the rcoal comparison

CXX      = g++
CXXFLAGS = -O2 -std=c++17 -Wall -Wextra -pthread
# Coalesce must not fuse multiplies and adds: its meshes would then differ from the Fortran's.
RCOAL_FLAGS = -ffp-contract=off

SRC = src/ses.cpp src/mesh.cpp src/clean.cpp src/io.cpp src/parallel.cpp
HDR = src/vec3.h src/ses.h src/mesh.h src/clean.h src/io.h src/parallel.h
RCOAL_SRC = src/coalesce.cpp
RCOAL_HDR = src/coalesce.h src/fortio.h

all: bin/sesmesh bin/rcoal bin/rcoal3

bin:
	mkdir -p bin

bin/sesmesh: src/sesmesh.cpp $(SRC) $(HDR) | bin
	$(CXX) $(CXXFLAGS) -o $@ src/sesmesh.cpp $(SRC)

bin/rcoal: src/rcoal.cpp $(RCOAL_SRC) $(RCOAL_HDR) | bin
	$(CXX) $(CXXFLAGS) $(RCOAL_FLAGS) -o $@ src/rcoal.cpp $(RCOAL_SRC)

# rcoal 3.0's messages by default (the program checks the name it was run as)
bin/rcoal3: bin/rcoal
	cp bin/rcoal $@

bin/test_sesmesh: test/test_sesmesh.cpp $(SRC) $(HDR) | bin
	$(CXX) $(CXXFLAGS) -Isrc -o $@ test/test_sesmesh.cpp $(SRC)

bin/test_coalesce: test/test_coalesce.cpp $(SRC) $(HDR) $(RCOAL_SRC) $(RCOAL_HDR) | bin
	$(CXX) $(CXXFLAGS) $(RCOAL_FLAGS) -Isrc -o $@ test/test_coalesce.cpp $(SRC) $(RCOAL_SRC)

test: bin/test_sesmesh bin/test_coalesce bin/rcoal bin/rcoal3
	bin/test_sesmesh
	bin/test_coalesce
	sh test/test_rcoal.sh

clean:
	rm -rf bin

.PHONY: all test clean
