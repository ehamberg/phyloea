OBJS=src/PhyloTree.o src/EvolutionModel.o

CXX=g++
CXXFLAGS=-Wall -pedantic -Werror -g
LIBS=
TESTLIBS=-lgtest -lpthread

TARGET=main
TESTTARGET=runtests

all: main tests

main: src/main.o $(OBJS)
	$(CXX) $(CXXFLAGS) src/main.o $(OBJS) -o $(TARGET) $(LIBS)

tests: tests/main.o $(OBJS)
	$(CXX) tests/main.o $(OBJS) -o $(TESTTARGET) $(TESTLIBS)

clean:
	rm -f */*.o
