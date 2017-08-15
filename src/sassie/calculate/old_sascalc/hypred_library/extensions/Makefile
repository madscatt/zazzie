CC=g++
CCFLAGS=-std=c++11 -c -fPIC -Wall -g
LDFLAGS=-std=c++11 -g

DIRSRC=src/
DIRINC=include/
DIRLIB=lib/
DIRTESTS=tests/cpp/

SRC=box.cc coor.cc pRDF.cc util.cc Atom.cc
LIB=prdf
SOURCES=$(addprefix $(DIRSRC),$(SRC))
HEADERS=$(addprefix $(DIRINC),$(SRC:.cc=.h))
INLINES=$(addprefix $(DIRINC),$(SRC:.cc=.icc))
OBJECTS=$(addprefix $(DIRLIB),$(SRC:.cc=.o))
LIBRARY=$(DIRLIB)lib$(LIB).a

#a.out: $(DIRTESTS)/main.cc $(LIBRARY)
#	$(CC) $(LDFLAGS) -L$(DIRLIB) -l$(LIB) -o $@ $<

$(LIBRARY): $(OBJECTS)
	$(AR) rv  $@ $^

$(OBJECTS): $(DIRLIB)%.o: $(DIRSRC)%.cc $(DIRINC)%.h $(DIRINC)%.icc
	$(CC) $(CCFLAGS) $< -o $@

.PHONY: clean
clean:
	rm -f $(OBJECTS)
	rm -f $(LIBRARY)
	rm -f a.out 
