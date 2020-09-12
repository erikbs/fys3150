# Vel clang om han finst
ifneq ($(shell command -v clang++ 2> /dev/null),)
  CXX := clang++
endif

CXXFLAGS += -std=c++17 -O2
CXXFLAGS += -I/opt/local/include
LDFLAGS += -larmadillo
LDFLAGS += -L/opt/local/lib

all: main

.PHONY: clean

main: main.cpp Tridiagonal.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o main main.cpp Tridiagonal.o

Tridiagonal.o: Tridiagonal.cpp Tridiagonal.hpp
	$(CXX) $(CXXFLAGS) -c -o Tridiagonal.o Tridiagonal.cpp

debug: CXXFLAGS += -DDEBUG
debug: all

clean:
	rm -rf *.o

distclean: clean
	rm -rf main genvec
