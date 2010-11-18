OBJS := $(patsubst %.cpp,%.o,$(wildcard src/*.cpp))
TESTOBJS := $(patsubst %.cpp,%.o,$(wildcard tests/*.cpp))

CXX=g++
CXXFLAGS=-Wall -pedantic -Isrc -g

LIBS=
TARGET=./main

TESTLIBS=-lgtest -lpthread
TESTTARGET=tests/runtests

EVTARGET=evaluator

all: main test

main: main.o $(OBJS)
	$(CXX) main.o $(OBJS) $(CXXFLAGS) -o $(TARGET) $(LIBS)

test: $(OBJS) $(TESTOBJS)
	$(CXX) $(TESTOBJS) $(OBJS) -o $(TESTTARGET) $(TESTLIBS)
	$(TESTTARGET) --gtest_shuffle 2> /dev/null

evaluator: evaluator.o $(OBJS)
	$(CXX) evaluator.o $(OBJS) $(CXXFLAGS) -o $(EVTARGET) $(LIBS)

memcheck:
	valgrind --tool=memcheck --leak-check=yes $(TARGET)

clean:
	rm -f */*.o *.o
