CPPC                 = gcc
CPPFLAGS             = -std=c++14 -O2 -Wall
LIBS = -liRRAM -lmpfr -lgmp -lpng
INCLUDES = -I./include

##########
# path to the library directory
LIB_DIR = ./lib

#########
# all cpp files
SRC_CPP_FILES := $(wildcard ./src/*/*.cpp)
SRC_OBJ_FILES := $(subst .cpp,.o,$(SRC_CPP_FILES))
INC_HPP_FILES := $(wildcard ./include/iRRAMx/*.hpp) $(wildcard ./include/iRRAMx/*/*.hpp)

##########
# target is the library file
MAIN = $(LIB_DIR)/libiRRAMx.a

all: $(MAIN)

$(MAIN): $(SRC_OBJ_FILES) lib $(INC_HPP_FILES)
	ar -r $(MAIN) $(SRC_OBJ_FILES)

##########
# create ./lib if it does not exist
lib:
	mkdir $(LIB_DIR)


##########
%.o: %.cpp
	$(CPPC) $(INCLUDES) $(CPPFLAGS)  -c $< -o $@


##########
# run test code

TESTDIR = ./test

test:
	g++ $(CPPFLAGS) $(INCLUDES) $(TESTDIR)/test.cpp -L./$(LIB_DIR) -liRRAMx $(LIBS) -o $(TESTDIR)/test



#########
# development purpose: create independent sub-library

##########
# paths to the source dirs of the sublibraries
LINEAR_DIR = ./src/linear
COMPACT_DIR = ./src/compact
POLYNOMIAL_DIR = ./src/polynomial
RANDOM_DIR = ./src/random
PLOT_DIR = ./src/plot

##########
# source and object files of each sublibrary
LINEAR_SRC_FILES := $(wildcard $(LINEAR_DIR)/*.cpp)
LINEAR_OBJ_FILES := $(patsubst $(LINEAR_DIR)/%.cpp,$(LINEAR_DIR)/%.o,$(LINEAR_SRC_FILES))

COMPACT_SRC_FILES := $(wildcard $(COMPACT_DIR)/*.cpp)
COMPACT_OBJ_FILES := $(patsubst $(COMPACT_DIR)/%.cpp,$(COMPACT_DIR)/%.o,$(COMPACT_SRC_FILES))

POLYNOMIAL_SRC_FILES := $(wildcard $(POLYNOMIAL_DIR)/*.cpp)
POLYNOMIAL_OBJ_FILES := $(patsubst $(POLYNOMIAL_DIR)/%.cpp,$(POLYNOMIAL_DIR)/%.o,$(POLYNOMIAL_SRC_FILES))

RANDOM_SRC_FILES := $(wildcard $(RANDOM_DIR)/*.cpp)
RANDOM_OBJ_FILES := $(patsubst $(RANDOM_DIR)/%.cpp,$(RANDOM_DIR)/%.o,$(RANDOM_SRC_FILES))

PLOT_SRC_FILES := $(wildcard $(PLOT_DIR)/*.cpp)
PLOT_OBJ_FILES := $(patsubst $(PLOT_DIR)/%.cpp,$(PLOT_DIR)/%.o,$(PLOT_SRC_FILES))

##########
# Sublibrary build for development
linear: $(LINEAR_OBJ_FILES) lib
	ar -r $(LIB_DIR)/libiRRAM_linear.a $(LINEAR_OBJ_FILES)

compact: $(COMPACT_OBJ_FILES) lib
	ar -r $(LIB_DIR)/libiRRAM_compact.a $(COMPACT_OBJ_FILES)

polynomial: $(POLYNOMIAL_OBJ_FILES) lib
	ar -r $(LIB_DIR)/libiRRAM_polynomial.a $(POLYNOMIAL_OBJ_FILES)

random: $(RANDOM_OBJ_FILES) lib
	ar -r $(LIB_DIR)/libiRRAM_random.a $(RANDOM_OBJ_FILES)

plot: $(PLOT_OBJ_FILES) lib
	ar -r $(LIB_DIR)/libiRRAM_plot.a $(PLOT_OBJ_FILES)

ALL_LIBS = $(wildcard $(LIB_DIR)/*.a)


##########
# compile test code that uses a sublibrary (for development use)
test_linear:
	g++ $(CPPFLAGS) $(INCLUDES) $(TESTDIR)/test_linear.cpp -L./$(LIB_DIR) -liRRAM_linear $(LIBS) -o $(TESTDIR)/test_linear

test_compact:
	g++ $(CPPFLAGS) $(INCLUDES) $(TESTDIR)/test_compact.cpp -L./$(LIB_DIR) -liRRAM_compact $(LIBS) -o $(TESTDIR)/test_compact

test_polynomial:
	g++ $(CPPFLAGS) $(INCLUDES) $(TESTDIR)/test_polynomial.cpp -L./$(LIB_DIR) -liRRAM_polynomial $(LIBS) -o $(TESTDIR)/test_polynomial

test_random:
	g++ $(CPPFLAGS) $(INCLUDES) $(TESTDIR)/test_random.cpp -L./$(LIB_DIR) -liRRAM_random $(LIBS) -o $(TESTDIR)/test_random

test_plot:
	g++ $(CPPFLAGS) $(INCLUDES) $(TESTDIR)/test_plot.cpp -L./$(LIB_DIR) -liRRAM_plot $(LIBS) -o $(TESTDIR)/test_plot


##########
# clean
.PHONY : clean test test_linear test_compact test_polynomial test_random test_plot
clean:
	-rm $(SRC_OBJ_FILES) $(ALL_LIBS)
