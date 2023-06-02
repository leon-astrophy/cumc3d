# Choose between CPU or GPU threading #
MULTICORE = -DGPU -acc #-DDEBUG -Minfo=accel
#MULTICORE = -mp #-DDEBUG -Minfo=accel

# Define F90 compiler #
F90=/home/leon/HDF5/bin/h5fc

# Define F90 flags #
F90FLAGS= -llapack -lblas -cpp $(MULTICORE) -fastsse -Mipa=fast -O4 -Mvect=simd -Mvect=levels:5 -Msmartalloc -Mstack_arrays -Mprefetch -cpp $(MULTICORE)
LDFLAGS= -llapack -lblas -cpp $(MULTICORE) -fastsse -Mipa=fast -O4 -Mvect=simd -Mvect=levels:5 -Msmartalloc -Mstack_arrays -Mprefetch -cpp $(MULTICORE)

# Define HDF5 library location
HDF5INCS=-I/home/leon/HDF5/include
HDF5LIBS=-L/home/leon/HDF5/lib
