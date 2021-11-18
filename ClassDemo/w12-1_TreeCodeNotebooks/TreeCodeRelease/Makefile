CXXFLAGS := -std=c++2a -march=native -fopenmp
OPT := -Ofast

DEFS =

LIBBASE := treelib
SHLIB := treelib.so

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
CXX := g++
LIBEXT := .so
SHFLAGS := -shared
endif

ifeq ($(UNAME_S),Darwin)
CXX = clang++
LIBEXT := .dylib
SHFLAGS := -dynamiclib
endif

SHLIB := $(addsuffix $(LIBEXT), $(LIBBASE))

SRC_DIR := src
LIB_DIR := lib

INC := -I$(SRC_DIR) -I$(LIB_DIR)

SRC := $(wildcard $(SRC_DIR)/*.cc)
LIB := $(wildcard $(LIB_DIR)/*.cc)

OBJ := $(SRC:.cc=.o)


# lib is not a file; it's only an easier-to-type target than treelib.so
.PHONY : lib
lib: $(SHLIB)

lib.o: $(SRC) $(LIB) lib.cc
	$(CXX) -o $@  $(DEFS) $(OPT) $(INC) $(CXXFLAGS) -fPIC -c lib.cc

$(SHLIB): lib.o
	$(CXX) $(SHFLAGS) -o $(SHLIB) $<

nbody: nbody.cc  $(SRC) $(LIB)
	$(CXX) -o $@  $(DEFS) $(OPT) $(INC) $(CXXFLAGS) $<

PHONY: clean
clean:
	rm -f *.o *~ a.out nbody treelib.so
