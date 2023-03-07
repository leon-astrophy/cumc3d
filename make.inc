F90=/home/leon/HDF5/bin/h5fc
F90FLAGS= -mp -fastsse -Mipa=fast -O4 -Mvect=simd -Mvect=levels:5 -Msmartalloc -Mstack_arrays -Mprefetch -Mconcur
LDFLAGS= -mp -fastsse -Mipa=fast -O4 -Mvect=simd -Mvect=levels:5 -Msmartalloc -Mstack_arrays -Mprefetch -Mconcur
HDF5INCS=-I/home/leon/HDF5/include
HDF5LIBS=-L/home/leon/HDF5/lib
