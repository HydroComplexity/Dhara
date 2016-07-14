# Compilers
MPICC=$(PREP) mpicc
MPILD=$(PREP) mpic++
NVCC=$(PREP) nvcc

# Source, Include, Object, and Library Directories
BINDIR=bin
SRCDIR=src
INCDIR=include
OBJDIR=src/obj

# Flags
CFLAGS      = -O3 -march=native -Wall
MPICFLAGS   = -I${MPI_HOME}/include
CUDACFLAGS  = -I${CUDA_PATH}/include
NETCDFFLAGS = -I${NETCDF4_HOME}/include
CUSPFLAGS   = -I ./include/cusplibrary-0.5.1/
EIGENFLAGS  = -I ./include

XCUDAFE		  := -Xcudafe "--diag_suppress=boolean_controlling_expr_is_constant"
GENCODE_SM13  := -gencode arch=compute_13,code=sm_13
GENCODE_SM20  := -gencode arch=compute_20,code=sm_20
GENCODE_SM30  := -gencode arch=compute_30,code=sm_30
GENCODE_SM35  := -gencode arch=compute_35,code=sm_35
GENCODE_FLAGS := $(GENCODE_SM20) $(GENCODE_SM30) $(GENCODE_SM35)

NVCCFLAGS     = -O3 $(GENCODE_FLAGS) -Xcompiler -march=native
CUDALDFLAGS   = -L${CUDA_PATH}/lib64 -lcudart
NETCDFLDFLAGS = -L${NETCDF4_HOME}/lib -lnetcdf

DHARA    = $(BINDIR)/dhara
BINARIES = $(DHARA)

# Look for all source files.
SRCS_CC = $(wildcard $(SRCDIR)/*.cc)
SRCS_CU = $(wildcard $(SRCDIR)/*.cu)

# Create an object for each source file and add object dir.
OBJECTS = $(addprefix $(OBJDIR)/, $(notdir $(SRCS_CC:.cc=.o) $(SRCS_CU:.cu=.cu.o)))

# Commands
all: build print $(BINARIES) done

$(OBJDIR)/%.cu.o: $(INCDIR)/main.h $(INCDIR)/cusplib.h $(SRCDIR)/%.cu Makefile
	$(NVCC) $(MPICFLAGS) $(NVCCFLAGS) $(XCUDAFE) $(CUSPFLAGS) $(EIGENFLAGS) $(NETCDFFLAGS) -c $(SRCDIR)/$*.cu -o $@

$(OBJDIR)/%.o: $(INCDIR)/*.h $(SRCDIR)/%.cc Makefile
	$(MPILD) $(MPICFLAGS) $(CFLAGS) $(CUDACFLAGS) $(EIGENFLAGS) $(NETCDFFLAGS) -c $(SRCDIR)/$*.cc -o $@

$(DHARA): $(OBJECTS) Makefile
	$(MPILD) $(CUDALDFLAGS) $(CUSPFLAGS) $(EIGENFLAGS) $(NETCDFLDFLAGS) -o $(DHARA) $(OBJECTS)

print:
	@echo ""
	@echo "Compling objects..."

# Check for the existence of sub directories.
# If not exist, then create these directories to dump data to.
# If exist do nothing
build:
	@echo ""
	@echo "Creating Directories..."
	if [ ! -d "$(OBJDIR)" ]; then mkdir $(OBJDIR);  fi
	if [ ! -d "$(BINDIR)" ]; then mkdir $(BINDIR);  fi

done: $(BINARIES)
	@echo ""
	@echo "Done!!!!"

clean:
	rm -rf $(OBJDIR)/*.o $(OBJDIR)/*.cu.o *~ $(BINARIES)
