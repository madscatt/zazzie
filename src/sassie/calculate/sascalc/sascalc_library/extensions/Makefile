# directories
DIRSRC=src/
DIRINC=include/
DIRLIB=lib/

ifndef _ARCH
    _ARCH := $(shell uname)
    export _ARCH
endif

INCDIR:=-I./include/
# compilers and flags
ifeq ($(_ARCH), Darwin)
    INCDIR+=-I/usr/local/cuda/include
    CC=g++
    NVCC=/usr/local/cuda/bin/nvcc
endif
ifeq ($(_ARCH), Linux)
    INCDIR+=-I/share/apps/local/cuda/include
    CC=/share/apps/local/bin/g++
    NVCC=/share/apps/local/cuda/bin/nvcc
endif

CCFLAGS=-fPIC -g -std=c++11

NVCCFLAGS=-Xcompiler -fPIC -DCUCC_ID="nVidia" -DBLD_PROCEDURE -DWITH_GSL -DHAVE_INLINE -DWITH_CUDA -DWITH_MPI -DPARALLEL -DWITH_PARALLEL_IO  $(CUDA_ARCH) -g -G -arch=sm_30

LIBNAME=sascalc
CUDA_LIBNAME=cudaSascalc

LIB=$(DIRLIB)lib$(LIBNAME).a
CUDA_LIB=$(DIRLIB)lib$(CUDA_LIBNAME).a

OBJS=SasCalc ScVars ScResults GV cudaGV
SOURCES=$(addsuffix .cpp,$(addprefix $(DIRSRC),$(OBJS)))
HEADERS=$(addsuffix .h,$(addprefix $(DIRINC),$(OBJS)))
OBJECTS=$(addsuffix .o,$(addprefix $(DIRLIB),$(OBJS)))

CUDA_OBJS=cudaKernel_GV wrapperCudaKernel
CUDA_SOURCES=$(addsuffix .cu,$(addprefix $(DIRSRC),$(CUDA_OBJS)))
CUDA_HEADERS=$(addsuffix .h,$(addprefix $(DIRINC),$(CUDA_OBJS)))
CUDA_OBJECTS=$(addsuffix .o,$(addprefix $(DIRLIB),$(CUDA_OBJS)))

# dependencies
all: $(LIB) $(CUDA_LIB) 

.PHONY: cpu
cpu: $(LIB)
$(LIB): $(OBJECTS)
	ar rv $@ $^

$(OBJECTS): $(DIRLIB)%.o: $(DIRSRC)%.cpp $(DIRINC)%.h
	$(CC) -c $< -o $@ $(CCFLAGS) $(INCDIR)

.PHONY: gpu
gpu: $(CUDA_LIB)
$(CUDA_LIB): $(CUDA_OBJECTS)
	ar rv $@ $^

$(CUDA_OBJECTS): $(DIRLIB)%.o: $(DIRSRC)%.cu $(DIRINC)%.h
	$(NVCC) -c $< -o $@ $(NVCCFLAGS) $(INCDIR)

.PHONY: clean
clean:
	rm -rf $(OBJECTS) $(CUDA_OBJECTS) $(LIB) $(CUDA_LIB)
