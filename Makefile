OBJS := $(patsubst %.cpp,%.o,$(wildcard src/*.cpp))
TESTOBJS := $(patsubst %.cpp,%.o,$(wildcard tests/*.cpp))

CXX=g++
CXXFLAGS=-Wall -pedantic -Isrc -g

LIBS=
TARGET=./master

TESTLIBS=-lgtest -lpthread
TESTTARGET=tests/runtests

EVTARGET=evaluator

.PHONY:all
all: master evaluator test

src/PhyloTree.o: src/PhyloTree.cpp src/PhyloTree.h
src/EASystem.o: src/EASystem.cpp src/EASystem.h
src/AMaxOperators.o: src/AMaxOperators.cpp src/AMaxOperators.h
src/EvolutionModel.o: src/EvolutionModel.cpp src/EvolutionModel.h
src/Fasta.o: src/Fasta.cpp src/Fasta.h
src/TreeOperators.o: src/TreeOperators.cpp src/TreeOperators.h
main.o: src/EAOperators.h src/EvolutionModel.h src/Fasta.h src/EASystem.h
mut.o: mut.cpp
recomb.o: recomb.cpp

.PHONY:master
master: main.o $(OBJS)
	@echo "$(OBJS)"
	$(CXX) main.o $(OBJS) $(CXXFLAGS) -o $(TARGET) $(LIBS)

.PHONY:mut
mut: mut.o $(OBJS)
	@echo "$(OBJS)"
	$(CXX) mut.o $(OBJS) $(CXXFLAGS) -o mut $(LIBS)

.PHONY:recomb
recomb: recomb.o $(OBJS)
	@echo "$(OBJS)"
	$(CXX) recomb.o $(OBJS) $(CXXFLAGS) -o recomb $(LIBS)

.PHONY:test
test: $(OBJS) $(TESTOBJS)
	$(CXX) $(TESTOBJS) $(OBJS) -o $(TESTTARGET) $(TESTLIBS)
	TERM=xterm $(TESTTARGET) --gtest_shuffle 2> /dev/null

.PHONY:evaluator
evaluator: evaluator.o $(OBJS)
	$(CXX) evaluator.o $(OBJS) $(CXXFLAGS) -o $(EVTARGET) $(LIBS)

.PHONY:memcheck
memcheck:
	valgrind --tool=memcheck --leak-check=yes $(TARGET)

.PHONY:docs
docs:
	doxygen Doxyfile


.PHONY:clean
clean:
	rm -f */*.o *.o core
