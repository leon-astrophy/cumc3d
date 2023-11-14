!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Written by Leung Shing Chi in 2016
! The subroutine setup the primitive and conservative equation indexes
! Notice that if you want your to add your own quantities, you need to specify them here
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SETUP_EQNS
USE DEFINITION
USE MHD_MODULE
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Number of equations !
INTEGER :: no_of_eq

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialize no_of_eq
no_of_eq = 0

! Initialize imax imin !
imin = 0
imax = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for buliding equation indexes !

WRITE(*,*) 'Now we arrange no_of_eq according'
WRITE(*,*) 'For each advectable scalar, no_of_eq increases by 1'
WRITE(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now we do the NM sector

! Set up minimum equation for NM !
imin = no_of_eq + 1

! NM density
no_of_eq = no_of_eq + 1
imax = no_of_eq
irho = no_of_eq
write(*,*) 'Make irho = ', no_of_eq

! NM vel-x
no_of_eq = no_of_eq + 1
imax = no_of_eq
ivx = no_of_eq
write(*,*) 'Make ivx = ', no_of_eq

! NM vel-y
no_of_eq = no_of_eq + 1
imax = no_of_eq
ivy = no_of_eq
write(*,*) 'Make ivy = ', no_of_eq

! NM vel-z
no_of_eq = no_of_eq + 1
imax = no_of_eq
ivz = no_of_eq
write(*,*) 'Make ivz = ', no_of_eq

! NM epsilon
no_of_eq = no_of_eq + 1
imax = no_of_eq
itau = no_of_eq
write(*,*) 'Make itau = ', no_of_eq

! Custom equations !
CALL CUSTOM_EQN

! Magnetic fields, Bx !
no_of_eq = no_of_eq + 1
imax = no_of_eq
ibx = no_of_eq
write(*,*) 'Make ibx = ', no_of_eq

! Magnetic fields, By !
no_of_eq = no_of_eq + 1
imax = no_of_eq
iby = no_of_eq
write(*,*) 'Make iby = ', no_of_eq

! Magnetic fields, Bz !
no_of_eq = no_of_eq + 1
imax = no_of_eq
ibz = no_of_eq
write(*,*) 'Make ibz = ', no_of_eq

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For electri field !
iex = imax + 1
iey = imax + 2
iez = imax + 3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Flip signs according to boundary conditions !

! Assume all variables are even !
bfac_xin(:) = 1
bfac_yin(:) = 1
bfac_zin(:) = 1
bfac_xout(:) = 1
bfac_yout(:) = 1
bfac_zout(:) = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Flip signs, x inner boundary !
IF(boundary_flag(1) == 2) THEN
	bfac_xin(ivx) = -1
	bfac_xin(ibx) = -1
ELSEIF(boundary_flag(1) == 3) THEN
	bfac_xin(ivx) = -1
	bfac_xin(ibx) = -1
	bfac_xin(ivz) = -1
	bfac_xin(ibz) = -1
ELSEIF(boundary_flag(1) == 4) THEN
	STOP 'Equatorial symmetry is not allowed for the x-boundary'
END IF

! Flip signs, x outer boundary !
IF(boundary_flag(2) == 2) THEN
	bfac_xout(ivx) = -1
	bfac_xout(ibx) = -1
ELSEIF(boundary_flag(2) == 3) THEN
	bfac_xout(ivx) = -1
	bfac_xout(ibx) = -1
	bfac_xout(ivz) = -1
	bfac_xout(ibz) = -1
ELSEIF(boundary_flag(2) == 4) THEN
	STOP 'Equatorial symmetry is not allowed for the x-boundary'
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Flip signs, y inner boundary !
IF(boundary_flag(3) == 2) THEN
	bfac_yin(ivy) = -1
	bfac_yin(iby) = -1
ELSEIF(boundary_flag(3) == 3) THEN
	bfac_yin(ivy) = -1
	bfac_yin(iby) = -1
	bfac_yin(ivz) = -1
	bfac_yin(ibz) = -1
ELSEIF(boundary_flag(3) == 4) THEN
	bfac_yin(ivy) = -1
	bfac_yin(ibx) = -1
	bfac_yin(ibz) = -1
END IF

! Flip signs, y outer boundary !
IF(boundary_flag(4) == 2) THEN
	bfac_yout(ivy) = -1
	bfac_yout(iby) = -1
ELSEIF(boundary_flag(4) == 3) THEN
	bfac_yout(ivy) = -1
	bfac_yout(iby) = -1
	bfac_yout(ivz) = -1
	bfac_yout(ibz) = -1
ELSEIF(boundary_flag(4) == 4) THEN
	bfac_yout(ivy) = -1
	bfac_yout(ibx) = -1
	bfac_yout(ibz) = -1
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Flip signs, z inner boundary !
IF(boundary_flag(5) == 2) THEN
	bfac_zin(ivz) = -1
	bfac_zin(ibz) = -1
ELSEIF(boundary_flag(5) == 3) THEN
	STOP 'Axial symmetry is not allowed for the z-boundary'
ELSEIF(boundary_flag(5) == 4) THEN
	bfac_zin(ivz) = -1
	bfac_zin(ibx) = -1
	bfac_zin(iby) = -1
END IF

! Flip signs, z outer boundary !
IF(boundary_flag(6) == 2) THEN
	bfac_zout(ivz) = -1
	bfac_zout(ibz) = -1
ELSEIF(boundary_flag(6) == 3) THEN
	STOP 'Axial symmetry is not allowed for the z-boundary'
ELSEIF(boundary_flag(6) == 4) THEN
	bfac_zout(ivz) = -1
	bfac_zout(ibx) = -1
	bfac_zout(iby) = -1
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For x-electric field !

! Ex at the x-inner boundary !
IF(bfac_xin(ivz)*bfac_xin(iby)*bfac_xin(ivy)*bfac_xin(ibz) < 0) THEN
	STOP 'Error in x-inner boundary flags for vz and by'
ELSE
  bfac_xin(iex) = bfac_xin(ivz)*bfac_xin(iby)
END IF

! Ex at the x-outer boundary !
IF(bfac_xout(ivz)*bfac_xout(iby)*bfac_xout(ivy)*bfac_xout(ibz) < 0) THEN
	STOP 'Error in x-outer boundary flags for vz and by'
ELSE
  bfac_xout(iex) = bfac_xout(ivz)*bfac_xout(iby)
END IF

! Ex at the y-inner boundary !
IF(bfac_yin(ivz)*bfac_yin(iby)*bfac_yin(ivy)*bfac_yin(ibz) < 0) THEN
	STOP 'Error in y-inner boundary flags for vz and by'
ELSE
  bfac_yin(iex) = bfac_yin(ivz)*bfac_yin(iby)
END IF

! Ex at the y-outer boundary !
IF(bfac_yout(ivz)*bfac_yout(iby)*bfac_yout(ivy)*bfac_yout(ibz) < 0) THEN
	STOP 'Error in y-outer boundary flags for vz and by'
ELSE
  bfac_yout(iex) = bfac_yout(ivz)*bfac_yout(iby)
END IF

! Ex at the z-inner boundary !
IF(bfac_zin(ivz)*bfac_zin(iby)*bfac_zin(ivy)*bfac_zin(ibz) < 0) THEN
	STOP 'Error in z-inner boundary flags for vz and by'
ELSE
  bfac_zin(iex) = bfac_zin(ivz)*bfac_zin(iby)
END IF

! Ex at the z-outer boundary !
IF(bfac_zout(ivz)*bfac_zout(iby)*bfac_zout(ivy)*bfac_zout(ibz) < 0) THEN
	STOP 'Error in z-outer boundary flags for vz and by'
ELSE
  bfac_zout(iex) = bfac_zout(ivz)*bfac_zout(iby)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For y-electric field !

! Ey at the x-inner boundary !
IF(bfac_xin(ivx)*bfac_xin(ibz)*bfac_xin(ivz)*bfac_xin(ibx) < 0) THEN
	STOP 'Error in x-inner boundary flags for vx and bz'
ELSE
  bfac_xin(iey) = bfac_xin(ivx)*bfac_xin(ibz)
END IF

! Ey at the x-outer boundary !
IF(bfac_xout(ivx)*bfac_xout(ibz)*bfac_xout(ivz)*bfac_xout(ibx) < 0) THEN
	STOP 'Error in x-outer boundary flags for vx and bz'
ELSE
  bfac_xout(iey) = bfac_xout(ivx)*bfac_xout(ibz)
END IF

! Ey at the y-inner boundary !
IF(bfac_yin(ivx)*bfac_yin(ibz)*bfac_yin(ivz)*bfac_yin(ibx) < 0) THEN
	STOP 'Error in y-inner boundary flags for vx and bz'
ELSE
  bfac_yin(iey) = bfac_yin(ivx)*bfac_yin(ibz)
END IF

! Ey at the y-outer boundary !
IF(bfac_yout(ivx)*bfac_yout(ibz)*bfac_yout(ivz)*bfac_yout(ibx) < 0) THEN
	STOP 'Error in y-outer boundary flags for vx and bz'
ELSE
  bfac_yout(iey) = bfac_yout(ivx)*bfac_yout(ibz)
END IF

! Ey at the z-inner boundary !
IF(bfac_zin(ivx)*bfac_zin(ibz)*bfac_zin(ivz)*bfac_zin(ibx) < 0) THEN
	STOP 'Error in z-inner boundary flags for vx and bz'
ELSE
  bfac_zin(iey) = bfac_zin(ivx)*bfac_zin(ibz)
END IF

! Ey at the z-outer boundary !
IF(bfac_zout(ivx)*bfac_zout(ibz)*bfac_zout(ivz)*bfac_zout(ibx) < 0) THEN
	STOP 'Error in z-outer boundary flags for vx and bz'
ELSE
  bfac_zout(iey) = bfac_zout(ivx)*bfac_zout(ibz)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For z-electric field !

! Ez at the x-inner boundary !
IF(bfac_xin(ivy)*bfac_xin(ibx)*bfac_xin(ivx)*bfac_xin(iby) < 0) THEN
	STOP 'Error in x-inner boundary flags for vy and bx'
ELSE
  bfac_xin(iez) = bfac_xin(ivy)*bfac_xin(ibx)
END IF

! Ez at the x-outer boundary !
IF(bfac_xout(ivy)*bfac_xout(ibx)*bfac_xout(ivx)*bfac_xout(iby) < 0) THEN
	STOP 'Error in x-outer boundary flags for vy and bx'
ELSE
  bfac_xout(iez) = bfac_xout(ivy)*bfac_xout(ibx)
END IF

! Ez at the y-inner boundary !
IF(bfac_yin(ivy)*bfac_yin(ibx)*bfac_yin(ivx)*bfac_yin(iby) < 0) THEN
	STOP 'Error in y-inner boundary flags for vy and bx'
ELSE
  bfac_yin(iez) = bfac_yin(ivy)*bfac_yin(ibx)
END IF

! Ez at the y-outer boundary !
IF(bfac_yout(ivy)*bfac_yout(ibx)*bfac_yout(ivx)*bfac_yout(iby) < 0) THEN
	STOP 'Error in y-outer boundary flags for vy and bx'
ELSE
  bfac_yout(iez) = bfac_yout(ivy)*bfac_yout(ibx)
END IF

! Ez at the z-inner boundary !
IF(bfac_zin(ivy)*bfac_zin(ibx)*bfac_zin(ivx)*bfac_zin(iby) < 0) THEN
	STOP 'Error in z-inner boundary flags for vy and bx'
ELSE
  bfac_zin(iez) = bfac_zin(ivy)*bfac_zin(ibx)
END IF

! Ez at the z-outer boundary !
IF(bfac_zout(ivy)*bfac_zout(ibx)*bfac_zout(ivx)*bfac_zout(iby) < 0) THEN
	STOP 'Error in z-outer boundary flags for vy and bx'
ELSE
  bfac_zout(iez) = bfac_zout(ivy)*bfac_zout(ibx)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for printing out !

WRITE(*,*)
WRITE(*,*) 'Finished setting up equation indexes'
WRITE(*,*) 'There are', no_of_eq, 'number of equations'
WRITE(*,*)

END SUBROUTINE SETUP_EQNS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine allocate arrays related to hydrodynamics variables !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BUILD_HYDRO
USE DEFINITION
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Core arrays for solving hyperbolic PDE

! Allocate arrays for NM !
ALLOCATE(l_rk(imin:imax,-2:nx+3,-2:ny+3,-2:nz+3))
ALLOCATE(prim(imin:imax,-2:nx+3,-2:ny+3,-2:nz+3))
ALLOCATE(cons(imin:imax,-2:nx+3,-2:ny+3,-2:nz+3))
ALLOCATE(u_old(imin:imax,-2:nx+3,-2:ny+3,-2:nz+3))

! Allocate flux arrays for NM !
ALLOCATE(sc(imin:imax,-2:nx+3,-2:ny+3,-2:nz+3))
ALLOCATE(flux(imin:imax,-2:nx+3,-2:ny+3,-2:nz+3))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For Grid variables !
ALLOCATE (x(-2:nx+3))
ALLOCATE (y(-2:ny+3))
ALLOCATE (z(-2:nz+3))
ALLOCATE (xF(-3:nx+3))
ALLOCATE (yF(-3:ny+3))
ALLOCATE (zF(-3:nz+3))
ALLOCATE (dx(-2:nx+3))
ALLOCATE (dy(-2:ny+3))
ALLOCATE (dz(-2:nz+3))
ALLOCATE (xbar(-2:nx+3))
ALLOCATE (dx_cb(-2:nx+3))
ALLOCATE (dsine(-2:ny+3))
ALLOCATE (dcose(-2:ny+3))
ALLOCATE (sine(-2:ny+3))
ALLOCATE (sinf(-3:ny+3))
ALLOCATE (vol(-2:nx+3,-2:ny+3,-2:nz+3))

! Internal energy !
ALLOCATE (epsilon(-2:nx+3,-2:ny+3,-2:nz+3))

! Hydrodynamic variables !
ALLOCATE (cs(-2:nx+3,-2:ny+3,-2:nz+3))
ALLOCATE (dpdrho(-2:nx+3,-2:ny+3,-2:nz+3))
ALLOCATE (dpdeps(-2:nx+3,-2:ny+3,-2:nz+3))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reconstruction weight for PPM !

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Build you custom variables !
CALL CUSTOM_HYDRO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE