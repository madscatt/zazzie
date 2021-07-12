# directories
DIRSRC=src/
DIRINC=include/
DIRLIB=lib/

INSTALLDIR = ../../
INSTALLDIRINC=../../include/
INSTALLDIRLIB=../../lib/

# mac / linux
#SYSINSTALLDIRINC=/usr/local/core_libraries/include/
#SYSINSTALLDIRLIB=/usr/local/core_libraries/lib/
#CC=g++

#cluster
SYSINSTALLDIRINC=/share/apps/local/core_libraries/include/
SYSINSTALLDIRLIB=/share/apps/local/core_libraries/lib/
CC=/share/apps/local/bin/g++

# compilers and flags
CCFLAGS=-fPIC -g -std=c++11

LIBNAME=overlap

LIB=$(DIRLIB)lib$(LIBNAME).a

OBJS=overlap
SOURCES=$(addsuffix .cpp,$(addprefix $(DIRSRC),$(OBJS)))
HEADERS=$(addsuffix .h,$(addprefix $(DIRINC),$(OBJS)))
HEADERS_INSTALL=$(addsuffix .h,$(addprefix $(INSTALLDIRINC),$(OBJS)))
HEADERS_SYSINSTALL=$(addsuffix .h,$(addprefix $(SYSINSTALLDIRINC),$(OBJS)))
OBJECTS=$(addsuffix .o,$(addprefix $(DIRLIB),$(OBJS)))
OBJECTS_INSTALL=$(addsuffix .o,$(addprefix $(INSTALLDIRLIB),$(OBJS)))
OBJECTS_SYSINSTALL=$(addsuffix .o,$(addprefix $(SYSINSTALLDIRLIB),$(OBJS)))

# dependencies
all: $(LIB)

$(LIB): $(OBJECTS)
	ar rv $@ $^

$(OBJECTS): $(DIRLIB)%.o: $(DIRSRC)%.cpp $(DIRINC)%.h
	$(CC) -c $< -o $@ $(CCFLAGS) -I$(DIRINC)

install: all
	cp -f $(DIRINC)*.h $(INSTALLDIRINC)
	cp -f $(DIRINC)*.h $(SYSINSTALLDIRINC)
	cp -f $(DIRLIB)*.o $(INSTALLDIRLIB)
	cp -f $(DIRLIB)*.a $(INSTALLDIRLIB)
	cp -f $(DIRLIB)*.o $(SYSINSTALLDIRLIB)
	cp -f $(DIRLIB)*.a $(SYSINSTALLDIRLIB)

.PHONY: clean
clean:
	rm -rf $(OBJECTS) $(LIB)
	rm -rf $(HEADERS_INSTALL)
	rm -rf $(HEADERS_SYSINSTALL)
	rm -rf $(OBJECTS_INSTALL)
	rm -rf $(OBJECTS_SYSINSTALL)

