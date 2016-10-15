DIRINC=include/
DIRSRC=src/
DIRLIB=lib/
DIRTESTS=test/

CC=g++
CCFLAGS=-std=c++11 -c -fPIC -w -g -I$(DIRINC) -I/usr/local/include
#CCFLAGS=-std=c++11 -c -fPIC -w -g -I$(DIRINC) -I/home/hailiang/work/tools/hdf5/include
LDFLAGS=-std=c++11 -g -w -I$(DIRINC)
#CCFLAGS=-c -fPIC -w -g -I$(DIRINC) -stdlib=libstdc++ -stdlib=libc++ -std=c++11
#LDFLAGS=-g -w -I$(DIRINC) -stdlib=libstdc++ -stdlib=libc++ -std=c++11

SRC=SasCalc.cpp ScVars.cpp ScResults.cpp GV.cpp
LIB=sascalc
HEADERS=$(addprefix $(DIRINC),$(SRC:.cpp=.h))
SOURCES=$(addprefix $(DIRSRC),$(SRC))
OBJECTS=$(addprefix $(DIRLIB),$(SRC:.cpp=.o))
LIBRARY=$(DIRLIB)lib$(LIB).a

#all: $(DIRTESTS)test.out $(DIRTESTS)test_pteros.out
all: $(LIBRARY)

$(DIRTESTS)test.out: $(DIRTESTS)test.cpp $(LIBRARY)
	$(CC) $< $(LDFLAGS) -L$(DIRLIB_BOOST) -L$(DIRLIB) -l$(LIB) -lboost_regex-mt -o $@

$(LIBRARY): $(OBJECTS)
	$(AR) rv  $@ $^

$(OBJECTS): $(DIRLIB)%.o: $(DIRSRC)%.cpp $(DIRINC)%.h
	$(CC) $(CCFLAGS) $< -o $@

.PHONY: clean
clean:
	rm -f $(OBJECTS)
	rm -f $(LIBRARY)
	rm -f $(DIRTESTS)test.out
