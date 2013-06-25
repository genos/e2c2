CXX=g++-4.7  # also works with clang++
CXXFLAGS=-c -Wall -Wextra -pedantic -O3 -m64 -std=c++11
LD=g++-4.7  # also works with clang++
LDFLAGS=-Wall -Wextra -lntl -lgmp
SOURCES=curves_test.cc points_test.cc key_demo.cc ecoh_echo.cc
EXECUTABLES=$(SOURCES:.cc=.out)

.PHONY: tidy clean

all: $(SOURCES) $(EXECUTABLES)

.SECONDARY:

%.out: %.o
	$(LD) $(LDFLAGS) $^ -o $@

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $^ -o $@

tidy:
	rm -f *.o

clean:
	rm -f *.out *.o
