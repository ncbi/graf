BASEDIR = ./
INCLUDES = -I. -I$(BASEDIR)

# set up compiler and options
CXX = g++
CXXFLAGS = -std=c++11 -pthread -g $(INCLUDES)

CXXTEST = ../cxxtest-4.3
CXXTESTFLAGS = -I$(CXXTEST)

HDIR = ./
SRCDIR = ./


#-----Suffix Rules---------------------------
# set up C++ suffixes and relationship between .C and .o files

.SUFFIXES: .cpp .C .cc

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<


#-----File Dependencies----------------------

SRC = Util.cpp XchrFpSnps.cpp BimFileSnps.cpp FamSubjects.cpp SexCheck.cpp SampleGenoDist.cpp FpSnps.cpp SnpInfo.cpp Relatives.cpp SampleGeno.cpp Pedigree.cpp GRAF.cpp

OBJ = $(addsuffix .o, $(basename $(SRC)))

graf: $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ)

# --------------- Unit Test Targets ----------------

TESTGEN = $(CXXTEST)/bin/cxxtestgen

TEST_SRC = Test_BimFileSnps.h Test_XchrFpSnps.h Test_FamSubjects.h

TEST_OBJ = $(addsuffix .o, $(basename $(TEST_SRC)))

graf_test:  graf_test.cpp $(OBJ)
	$(CXX) $(CXXFLAGS) $(CXXTESTFLAGS) -o $@ $(filter-out GRAF.o, $(OBJ)) graf_test.cpp

graf_test.cpp: $(TEST_SRC)
	$(TESTGEN) --error-printer -o $@ $(TEST_SRC)

test: graf_test
	./graf_test

#-----Other stuff----------------------------

# makedepend does not work anymore
# depend:
# 	makedepend -- $(CXXFLAGS) -- -Y $(SRC) $(TEST_SRC)

clean:
	rm -f $(OBJ) $(TEST_OBJ) *~
	rm -f graf_test.cpp

# DO NOT DELETE

Util.o: Graf.h
XchrFpSnps.o: Graf.h XchrFpSnps.h
BimFileSnps.o: Graf.h XchrFpSnps.h BimFileSnps.h
FamSubjects.o: Graf.h FamSubjects.h
SexCheck.o: Graf.h FamSubjects.h SexCheck.h XchrFpSnps.h BimFileSnps.h
SampleGenoDist.o: Graf.h
FpSnps.o: Graf.h
SnpInfo.o: Graf.h
Relatives.o: Graf.h
SampleGeno.o: Graf.h
Pedigree.o: Graf.h
GRAF.o: Graf.h XchrFpSnps.h BimFileSnps.h FamSubjects.h SexCheck.h
Test_BimFileSnps.o: Graf.h XchrFpSnps.h BimFileSnps.h
Test_XchrFpSnps.o: Graf.h XchrFpSnps.h
Test_FamSubjects.o: Graf.h FamSubjects.h
