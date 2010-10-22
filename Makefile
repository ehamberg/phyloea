OBJS=src/main.o src/PhyloTree.o src/EvolutionModel.o

CXX=g++
CXXFLAGS=-Wall -pedantic -Werror -g
LIBS=

TARGET=test

all: $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET) $(LIBS)

clean:
	rm -f src/*.o
