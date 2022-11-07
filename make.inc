F90=nvfortran
F90FLAGS= -fastsse -Mipa=fast,inline -O4 -Mvect=simd -Msmartalloc -Mstack_arrays -Mprefetch -Mconcur -tp=px -stdpar=gpu,multicore -llapack -lblas -mp
MODINC="-I ./"
LDFLAGS= -fastsse -Mipa=fast,inline -O4 -Mvect=simd -Msmartalloc -Mstack_arrays -Mprefetch -Mconcur -tp=px -stdpar=gpu,multicore -llapack -lblas -mp
HDF5INCS=-I/home/leon/HDF5/include
HDF5LIBS=-L/home/leon/HDF5/lib -lhdf5 -lhdf5_fortran -lhdf5  

#must make clean then make if changed these pre-processor variables
HAVE_NUC_EOS=1
