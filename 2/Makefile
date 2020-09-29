# Vel clang om han finst
ifneq ($(shell command -v clang++ 2> /dev/null),)
  CXX := clang++
endif

CXXFLAGS += -std=c++17 
#-O2
CXXFLAGS += -I/opt/local/include
LDFLAGS += -larmadillo
LDFLAGS += -L/opt/local/lib

all: main test

.PHONY: clean

test: test.cpp Jacobi.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -DCATCH_CONFIG_DISABLE_EXCEPTIONS -o test test.cpp Jacobi.o

main: main.cpp Jacobi.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o main main.cpp Jacobi.o

Jacobi.o: Jacobi.cpp Jacobi.hpp
	$(CXX) $(CXXFLAGS) -c -o Jacobi.o Jacobi.cpp

plots: plot.py
	chmod a+x plot.py
	./plot.py

debug: CXXFLAGS += -DDEBUG
debug: all

clean:
	rm -rf *.o

distclean: clean
	rm -rf main test
