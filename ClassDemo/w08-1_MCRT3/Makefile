#
# Makefile for mcrt Monte Carlo transport code
#
# "make lib" to make the shared object library for python
# "make exec" to make executable using the C++ main in mcrt3.cc
#

UNAME = $(shell uname -s)
$(info OS name is ${UNAME})
$(info default C++ compiler $$CXX is ${CXX})

ifeq ($(UNAME),Linux)
	SHFLAGS = -shared
	SHLIB = libmcrt.so
	#CXX = g++
else
ifeq ($(UNAME),Monterey)
	SHFLAGS = -dynamiclib
	SHLIB = libmcrt.dylib
	#CXX = clang++
else
$(error Cannot determine OS)
endif
endif

# "phony" dependency as a convenient name
.PHONY: exec
exec: mcrt

# "phony" dependency as a convenient name
.PHONY: lib
lib: $(SHLIB)

# make the executable if any of these files are more recent than mcrt
mcrt: mcrt3.cc random.cc tallys.cc misc.cc STimer.cc vformat.cc
	$(CXX) -O3 -fopenmp mcrt3.cc -o mcrt

# make the shared library
$(SHLIB): mcrt3.o
	$(CXX) $(SHFLAGS)  -o $(SHLIB) mcrt3.o

# compile for the shared library
# defining PYTHONLIB excludes main in mcrt3.cc
mcrt3.o: mcrt3.cc random.cc tallys.cc misc.cc STimer.cc vformat.cc
	$(CXX) -std=c++2a -O3 -DPYTHONLIB -c -fPIC -fopenmp -o mcrt3.o mcrt3.cc

.PHONY: clean
clean:
	rm -f mcrt3.o $(SHLIB) mcrt
