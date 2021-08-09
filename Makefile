CXX=g++
CXXFLAGS=-std=c++17 -Wall -Wextra -O2

ising: main.o ising.o
	$(CXX) $(CXXFLAGS) main.o ising.o -o ising

main.o: main.cpp ising.hpp
	$(CXX) $(CXXFLAGS) -c main.cpp

ising.o: ising.cpp ising.hpp
	$(CXX) $(CXXFLAGS) -c ising.cpp

clean:
	rm *.o

.PHONY: clean
