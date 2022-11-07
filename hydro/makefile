include ../make.inc

SOURCES=definition.f90 \
	rsolver_module.f90 \
	weno_module.f90 \
	ppm_module.f90 \
	main.f90 \
	alphasplit.f90 \
	boundary.f90 \
	checkrho.f90 \
	dataio.f90 \
	findpressure.f90 \
	fromrvetou.f90 \
	getepsilon.f90 \
	getgrid.f90 \
	initial.f90 \
	riemann1d.f90 \
	riemann2d.f90 \
	riemann3d.f90 \
	rungekutta.f90 \
	setup.f90 \
	soundspeed.f90 \
	spatial.f90 \
	update.f90 

OBJECTS=$(SOURCES:.f90=.o )

EXTRAINCS += $(HDF5INCS)
EXTRAOBJECTS += $(HDF5LIBS)

CUMC3D: definition.o $(OBJECTS) $(EXTRADEPS)  
	$(F90) $(LDFLAGS) -o CUMC3D $(run) $(OBJECTS) $(EXTRAOBJECTS)

$(OBJECTS): %.o: %.f90 
	$(F90) $(F90FLAGS) $(DEFS) $(EXTRAINCS) $(MODINC) -c $< -o $@

definition.mod: definition.f90
	$(F90) -C definition.f90

rsolver.mod: rsolver_module.f90
	$(F90) -C rsolver_module.f90

weno.mod: weno_module.f90
	$(F90) -C weno_module.f90

clean:
	rm -rf CUMC3D
	rm -rf *.o 	
	rm -rf *.mod
	rm tmp.txt
cleanfile:
	rm -rf ../plot/*.png
	rm -rf ../outfile/*.hdf5
