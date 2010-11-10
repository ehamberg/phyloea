OBJS := $(subst src/main.o,, $(patsubst %.cpp,%.o,$(wildcard src/*.cpp)))
TESTOBJS := $(patsubst %.cpp,%.o,$(wildcard tests/*.cpp))

CXX=g++
CXXFLAGS=-Wall -pedantic -Isrc -g

LIBS=
TARGET=./main

TESTLIBS=-lgtest -lpthread
TESTTARGET=tests/runtests

all: src/main.o $(OBJS) $(TESTOBJS)
	$(CXX) src/main.o $(OBJS) $(CXXFLAGS) -o $(TARGET) $(LIBS)
	$(CXX) $(TESTOBJS) $(OBJS) -o $(TESTTARGET) $(TESTLIBS)
	$(TESTTARGET) --gtest_shuffle 2> /dev/null

memcheck:
	valgrind --tool=memcheck --leak-check=yes $(TARGET)

clean:
	rm -f */*.o
