F90=nvfortran
F90FLAGS= -mp -fastsse -Mipa=fast -O4 -Mvect=simd -Mvect=levels:5 -Msmartalloc -Mstack_arrays -Mprefetch -Mconcur
LDFLAGS= -mp -fastsse -Mipa=fast -O4 -Mvect=simd -Mvect=levels:5 -Msmartalloc -Mstack_arrays -Mprefetch -Mconcur
HDF5INCS=-I/home/leon/HDF5/include
HDF5LIBS=-L/home/leon/HDF5/lib -lhdf5 -lhdf5_fortran
