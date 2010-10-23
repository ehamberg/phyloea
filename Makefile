OBJS := $(subst src/main.o,, $(patsubst %.cpp,%.o,$(wildcard src/*.cpp)))
TESTOBJS := $(patsubst %.cpp,%.o,$(wildcard tests/*.cpp))

CXX=g++
CXXFLAGS=-Wall -pedantic -Werror -Isrc -g

LIBS=
TARGET=main

TESTLIBS=-lgtest -lpthread
TESTTARGET=runtests

all: src/main.o $(OBJS) $(TESTOBJS)
	$(CXX) src/main.o $(OBJS) $(CXXFLAGS) -o $(TARGET) $(LIBS)
	$(CXX) $(TESTOBJS) $(OBJS) -o $(TESTTARGET) $(TESTLIBS)
	./runtests --gtest_shuffle

clean:
	rm -f */*.o
