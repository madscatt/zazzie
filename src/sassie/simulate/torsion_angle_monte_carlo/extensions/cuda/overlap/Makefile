# directories
DIRSRC=src/
DIRINC=include/
DIRLIB=lib/


INSTALLDIRINC=../../include/
INSTALLDIRLIB=../../lib/

#mac / linux

SYSINSTALLDIRINC=/usr/local/core_libraries/include/
SYSINSTALLDIRLIB=/usr/local/core_libraries/lib/
INCDIR=-I/usr/local/cuda-12.1/include

#cluster
#SYSINSTALLDIRINC=/share/apps/local/core_libraries/include/
#SYSINSTALLDIRLIB=/share/apps/local/core_libraries/lib/
#INCDIR=-I/share/apps/local/cuda-6.5/include

# compilers and flags
CC=g++
CCFLAGS=-fPIC -g -std=c++11

NVCC=nvcc
NVCCFLAGS=-Xcompiler -fPIC -DCUCC_ID="nVidia" -DBLD_PROCEDURE -DWITH_GSL -DHAVE_INLINE -DWITH_CUDA -DWITH_MPI -DPARALLEL -DWITH_PARALLEL_IO  $(CUDA_ARCH) -g -G -arch=sm_30

CUDA_LIBNAME=cudaOverlap

CUDA_LIB=$(DIRLIB)lib$(CUDA_LIBNAME).a

CUDA_OBJS=cudaKernel_overlap wrapperCudaKernel
CUDA_SOURCES=$(addsuffix .cu,$(addprefix $(DIRSRC),$(CUDA_OBJS)))
CUDA_HEADERS=$(addsuffix .h,$(addprefix $(DIRINC),$(CUDA_OBJS)))
CUDA_HEADERS_INSTALL=$(addsuffix .h,$(addprefix $(INSTALLDIRINC),$(CUDA_OBJS)))
CUDA_HEADERS_SYSINSTALL=$(addsuffix .h,$(addprefix $(SYSINSTALLDIRINC),$(CUDA_OBJS)))
CUDA_OBJECTS=$(addsuffix .o,$(addprefix $(DIRLIB),$(CUDA_OBJS)))
CUDA_OBJECTS_INSTALL=$(addsuffix .o,$(addprefix $(INSTALLDIRLIB),$(CUDA_OBJS)))
CUDA_OBJECTS_SYSINSTALL=$(addsuffix .o,$(addprefix $(SYSINSTALLDIRLIB),$(CUDA_OBJS)))

# dependencies
all: $(CUDA_LIB)

$(CUDA_LIB): $(CUDA_OBJECTS)
	ar rv $@ $^

$(CUDA_OBJECTS): $(DIRLIB)%.o: $(DIRSRC)%.cu $(DIRINC)%.h
	$(NVCC) -c $< -o $@ $(NVCCFLAGS) -I$(DIRINC) $(INCDIR)

install: all
	cp -f $(DIRINC)*.h $(INSTALLDIRINC)
	cp -f $(DIRINC)*.h $(SYSINSTALLDIRINC)
	cp -f $(DIRLIB)*.o $(INSTALLDIRLIB)
	cp -f $(DIRLIB)*.o $(SYSINSTALLDIRLIB)
	cp -f $(DIRLIB)*.a $(INSTALLDIRLIB)
	cp -f $(DIRLIB)*.a $(SYSINSTALLDIRLIB)

.PHONY: clean
clean:
	rm -rf $(CUDA_OBJECTS) $(CUDA_LIB)
	rm -rf $(CUDA_HEADERS_INSTALL)
	rm -rf $(CUDA_HEADERS_SYSINSTALL)
	rm -rf $(CUDA_OBJECTS_INSTALL)
	rm -rf $(CUDA_OBJECTS_SYSINSTALL)
