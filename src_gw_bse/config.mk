# Define the following paths in order to build gw_bse
SRCDIR = /Users/kayahan/Documents/GitHub/OpenAtom_main/OpenAtom/src_gw_bse
BUILDDIR = $(SRCDIR)/build
# CHARMDIR = /Users/kayahan/Software/charm-aa/mpi-darwin-x86_64
CHARMDIR = /Users/kayahan/Software/charm-v7.0.0/mpi-darwin-x86_64-mpicxx

CHARMC = $(CHARMDIR)/bin/charmc

# Compile and link flags
CFLAGS = -g -O3 
LFLAGS = -module CkMulticast -mpi -nomain 

# Path to the build of the fftw3 library
FFTW3 = /opt/local
CFLAGS += -I$(FFTW3)/include
LFLAGS += -L$(FFTW3)/lib -lfftw3 -lm -lz

# Path to LAPACK if using LAPACK
# FOR MACS
LAPACKLIB = -L/opt/local/lib -lopenblas -lgomp -lgfortran
LAPACKINC = -I/opt/local/include/lapack -I/opt/local/include -lgfortran

SCALAPACKLIB = /Users/kayahan/Software/scalapack-2.1.0/libscalapack.a


# FOR LAB MACHINES
# LAPACKLIB = -llapack -lblas
# LAPACKINC =
# FOR ANL MACHINES
# LAPACKLIB = -L/soft/libraries/alcf/current/xl/CBLAS/lib -L/soft/libraries/alcf/current/xl/BLAS/lib -L/soft/compilers/ibmcmp-nov2012/xlf/bg/14.1/lib64 -lcblas -lblas -lxlf90_r -lxlopt -lxl -lxlfmath
# LAPACKINC = -I/soft/libraries/alcf/current/xl/CBLAS/include/

HDF5LIB =
#-L/usr/local/Cellar/hdf5/1.12.1/lib -lhdf5 -lhdf5_cpp
HDF5INC =
#-I/usr/local/Cellar/hdf5/1.12.1/include
CFLAGS += -DUSE_LAPACK -DUSE_FORTRAN_UNDERSCORE -DUSE_ZGEMM $(LAPACKINC) $(HDF5INC) 
LFLAGS += $(LAPACKLIB) $(SCALAPACKLIB) $(HDF5LIB)

# Flags for enabling CkLoop
CFLAGS += -DUSE_CKLOOP
LFLAGS += -module CkLoop

# The sub directories used for compilation
DIRS=configuration main fft matrix matmul states diagonalizer utils 
INCLUDES=$(addprefix -I$(SRCDIR)/, $(DIRS))
CFLAGS+=$(INCLUDES)


