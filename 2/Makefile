# Vel clang om han finst
ifneq ($(shell command -v clang++ 2> /dev/null),)
  CXX := clang++
endif

CXXFLAGS += -std=c++17 -O2
CXXFLAGS += -I/opt/local/include
LDFLAGS += -larmadillo
LDFLAGS += -L/opt/local/lib

all: jacobi test schrod

.PHONY: clean

test: test.cpp Jacobi.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -DCATCH_CONFIG_DISABLE_EXCEPTIONS -o test test.cpp Jacobi.o

jacobi: main.cpp Jacobi.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o jacobi main.cpp Jacobi.o

schrod: schrod.cpp Jacobi.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o schrod schrod.cpp Jacobi.o

schrod2: schrod2.cpp Jacobi.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o schrod2 schrod2.cpp Jacobi.o

Jacobi.o: Jacobi.cpp Jacobi.hpp
	$(CXX) $(CXXFLAGS) -c -o Jacobi.o Jacobi.cpp

plots: plot.py
	chmod a+x plot.py schrod.py
	./plot.py
	./schrod.py

debug: CXXFLAGS += -DDEBUG
debug: all

clean:
	rm -rf *.o

distclean: clean
	rm -rf jacobi schrod schrod2 test
