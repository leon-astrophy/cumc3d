include ../make.inc

SOURCES=definition.f90 \
	main.f90 \
	initial.f90 \
	setup.f90

OBJECTS=$(SOURCES:.f90=.o )

EXTRAINCS += $(HDF5INCS)
EXTRAOBJECTS += $(HDF5LIBS)

CUMC3D: definition.o $(OBJECTS) $(EXTRADEPS)  
	$(F90) $(LDFLAGS) -o CUMC3D $(run) $(OBJECTS) $(EXTRAOBJECTS)

$(OBJECTS): %.o: %.f90 
	$(F90) $(F90FLAGS) $(DEFS) $(EXTRAINCS) $(MODINC) -c $< -o $@

definition.mod: definition.f90
	$(F90) -C definition.f90

clean:
	rm -rf CUMC3D
	rm -rf *.o 	
	rm -rf *.mod
	rm tmp.txt