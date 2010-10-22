OBJS := $(subst src/main.o,, $(patsubst %.cpp,%.o,$(wildcard src/*.cpp)))
TESTOBJS := $(patsubst %.cpp,%.o,$(wildcard tests/*.cpp))

CXX=g++
CXXFLAGS=-Wall -pedantic -Werror -Isrc -g

LIBS=
TARGET=main

TESTLIBS=-lgtest -lpthread
TESTTARGET=runtests

all: main tests

main: src/main.o $(OBJS)
	$(CXX) src/main.o $(OBJS) $(CXXFLAGS) -o $(TARGET) $(LIBS)

tests: $(TESTOBJS) $(OBJS)
	$(CXX) $(TESTOBJS) $(OBJS) -o $(TESTTARGET) $(TESTLIBS)

clean:
	rm -f */*.o
