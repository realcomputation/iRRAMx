CPPC                 = gcc
STANDARD_FLAGS       = -std=c++14
# WARN_AS_ERRORS_FLAGS = -pedantic-errors\
#                        -Wall\
#                        -Wextra\
#                        -Werror\
#                        -Wconversion
# DEBUG_FLAGS          = -g3
# OPT_FLAGS            = -0O
# NO_LINKER_FLAGS      = -c

CPPFLAGS             = -std=c++14 -O2 -Wall 

# CXX = g++
# # CXXCPP = g++ -E -std=c++11
# CPPFLAGS = -I./include
# CXXFLAGS = -g -O2 -std=c++14 -mpc80
# LDFLAGS = -Xlinker -rpath -Xlinker ./lib
# LDLIBS = -L./lib -lmpfr -lgmp -lm -lpthread -lirram -lpng 


INCLUDES = -I./include

# SRCS     = ./src/compact/compact.cpp ./src/compact/compact2.cpp

OBJS     = ./src/compact/compact.o

LIBS = -lmpfr -lgmp -lm -lirram
MAIN = libtest.a # static library


# linear:
# 	gcc -c linear.cpp -liRRAM

all: $(MAIN)
	@echo $(MAIN) has been compiled!

$(MAIN): $(OBJS)
	ar -r $(MAIN) $(OBJS) 

%.o: %.cpp
	$(CPPC) $(INCLUDES) $(CPPFLAGS) $(LIBS) -c $< -o $@

test:
	gcc -I./include -std=c++14 -O2 -Wall -lmpfr -lgmp -lm -lirram test.cpp -o test
	# $(CPPC) $(INCLUDES) $(CPPFLAGS) $(LIBS) -o test test.cpp
# depend: $(SRCS)
# 	makedepend $(INCLUDES) $^