LD = $(CXX)
CXXFLAGS = -Wall -Wextra -pedantic -O3 -std=c++11
LDFLAGS = -lntl -lgmp
SOURCES = curves_test.cc points_test.cc key_demo.cc ecoh_echo.cc
EXECUTABLES = $(SOURCES:.cc=.out)

.PHONY: tidy clean

all: $(SOURCES) $(EXECUTABLES)

test: $(SOURCES) $(EXECUTABLES)
	./curves_test.out
	./points_test.out
	./key_demo.out
	./ecoh_echo.out 1234 secretPASSWORD

.SECONDARY:

%.out: %.o
	$(LD) $(LDFLAGS) $^ -o $@

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $^ -o $@

tidy:
	rm -f *.o

clean:
	rm -f *.out *.o
