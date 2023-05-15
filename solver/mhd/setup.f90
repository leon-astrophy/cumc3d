!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Written by Leung Shing Chi in 2016
! The subroutine setup the primitive and conservative equation indexes
! Notice that if you want your to add your own quantities, you need to specify them here
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SETUP_EQN
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
imin1 = 0
imax1 = 0
imin2 = 0
imax2 = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for buliding equation indexes !

WRITE(*,*) 'Now we arrange no_of_eq according'
WRITE(*,*) 'For each advectable scalar, no_of_eq increases by 1'
WRITE(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now we do the NM sector

! Set up minimum equation for NM !
imin2 = no_of_eq + 1

! NM density
no_of_eq = no_of_eq + 1
imax2 = no_of_eq
irho2 = no_of_eq
write(*,*) 'Make irho2 = ', no_of_eq

! NM vel-x
no_of_eq = no_of_eq + 1
imax2 = no_of_eq
ivel2_x = no_of_eq
write(*,*) 'Make ivel2_x = ', no_of_eq

! NM vel-y
no_of_eq = no_of_eq + 1
imax2 = no_of_eq
ivel2_y = no_of_eq
write(*,*) 'Make ivel2_y = ', no_of_eq

! NM vel-z
no_of_eq = no_of_eq + 1
imax2 = no_of_eq
ivel2_z = no_of_eq
write(*,*) 'Make ivel2_z = ', no_of_eq

! NM epsilon
no_of_eq = no_of_eq + 1
imax2 = no_of_eq
itau2 = no_of_eq
write(*,*) 'Make itau2 = ', no_of_eq

! Dual energies !
IF(dual_energy) THEN
	no_of_eq = no_of_eq + 1
	imax2 = no_of_eq
	ieps2 = no_of_eq
	write(*,*) 'Make ieps2 = ', no_of_eq
END IF

! Custom equations !
CALL CUSTOM_EQN

! Magnetic fields, Bx !
no_of_eq = no_of_eq + 1
imax2 = no_of_eq
ibx = no_of_eq
write(*,*) 'Make ibx = ', no_of_eq

! Magnetic fields, By !
no_of_eq = no_of_eq + 1
imax2 = no_of_eq
iby = no_of_eq
write(*,*) 'Make iby = ', no_of_eq

! Magnetic fields, Bz !
no_of_eq = no_of_eq + 1
imax2 = no_of_eq
ibz = no_of_eq
write(*,*) 'Make ibz = ', no_of_eq

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
	bfac_xin(ivel2_x) = -1
	bfac_xin(ibx) = -1
ELSEIF(boundary_flag(1) == 3) THEN
	bfac_xin(ivel2_x) = -1
	bfac_xin(ibx) = -1
	bfac_xin(ivel2_z) = -1
	bfac_xin(ibz) = -1
ELSEIF(boundary_flag(1) == 4) THEN
	STOP 'Equatorial symmetry is not allowed for the x-boundary'
END IF

! Flip signs, x outer boundary !
IF(boundary_flag(2) == 2) THEN
	bfac_xout(ivel2_x) = -1
	bfac_xout(ibx) = -1
ELSEIF(boundary_flag(2) == 3) THEN
	bfac_xout(ivel2_x) = -1
	bfac_xout(ibx) = -1
	bfac_xout(ivel2_z) = -1
	bfac_xout(ibz) = -1
ELSEIF(boundary_flag(2) == 4) THEN
	STOP 'Equatorial symmetry is not allowed for the x-boundary'
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Flip signs, y inner boundary !
IF(boundary_flag(3) == 2) THEN
	bfac_yin(ivel2_y) = -1
	bfac_yin(iby) = -1
ELSEIF(boundary_flag(3) == 3) THEN
	bfac_yin(ivel2_y) = -1
	bfac_yin(iby) = -1
	bfac_yin(ivel2_z) = -1
	bfac_yin(ibz) = -1
ELSEIF(boundary_flag(3) == 4) THEN
	bfac_yin(ivel2_y) = -1
	bfac_yin(ibx) = -1
	bfac_yin(ibz) = -1
END IF

! Flip signs, y outer boundary !
IF(boundary_flag(4) == 2) THEN
	bfac_yout(ivel2_y) = -1
	bfac_yout(iby) = -1
ELSEIF(boundary_flag(4) == 3) THEN
	bfac_yout(ivel2_y) = -1
	bfac_yout(iby) = -1
	bfac_yout(ivel2_z) = -1
	bfac_yout(ibz) = -1
ELSEIF(boundary_flag(4) == 4) THEN
	bfac_yout(ivel2_y) = -1
	bfac_yout(ibx) = -1
	bfac_yout(ibz) = -1
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Flip signs, z inner boundary !
IF(boundary_flag(5) == 2) THEN
	bfac_zin(ivel2_z) = -1
	bfac_zin(ibz) = -1
ELSEIF(boundary_flag(5) == 3) THEN
	STOP 'Axial symmetry is not allowed for the z-boundary'
ELSEIF(boundary_flag(5) == 4) THEN
	bfac_zin(ivel2_z) = -1
	bfac_zin(ibx) = -1
	bfac_zin(iby) = -1
END IF

! Flip signs, z outer boundary !
IF(boundary_flag(6) == 2) THEN
	bfac_zout(ivel2_z) = -1
	bfac_zout(ibz) = -1
ELSEIF(boundary_flag(6) == 3) THEN
	STOP 'Axial symmetry is not allowed for the z-boundary'
ELSEIF(boundary_flag(6) == 4) THEN
	bfac_zout(ivel2_z) = -1
	bfac_zout(ibx) = -1
	bfac_zout(iby) = -1
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Flip signs according to boundary conditions, for electric fields !

! initialize !
efac_xin(:) = 1
efac_xout(:) = 1
efac_yin(:) = 1
efac_yout(:) = 1
efac_zin(:) = 1
efac_zout(:) = 1

! initialize !
ecell_xin(:) = 1
ecell_xout(:) = 1
ecell_yin(:) = 1
ecell_yout(:) = 1
ecell_zin(:) = 1
ecell_zout(:) = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Flip signs, x inner boundary !
IF(boundary_flag(1) == 2) THEN
	efac_xin(iyx) = -1
	efac_xin(iyz) = -1
	efac_xin(izx) = -1
	efac_xin(izy) = -1
	ecell_xin(iy) = -1
	ecell_xin(iz) = -1
ELSEIF(boundary_flag(1) == 3) THEN
	efac_xin(ixy) = -1
	efac_xin(ixz) = -1
	efac_xin(izx) = -1
	efac_xin(izy) = -1
	ecell_xin(ix) = -1
	ecell_xin(iz) = -1
ELSEIF(boundary_flag(1) == 4) THEN
	STOP 'Equatorial symmetry is not allowed for the x-boundary'
END IF

! Flip signs, x outer boundary !
IF(boundary_flag(2) == 2) THEN
	efac_xout(iyx) = -1
	efac_xout(iyz) = -1
	efac_xout(izx) = -1
	efac_xout(izy) = -1
	ecell_xout(iy) = -1
	ecell_xout(iz) = -1
ELSEIF(boundary_flag(2) == 3) THEN
	efac_xout(ixy) = -1
	efac_xout(ixz) = -1
	efac_xout(izx) = -1
	efac_xout(izy) = -1
	ecell_xout(ix) = -1
	ecell_xout(iz) = -1
ELSEIF(boundary_flag(2) == 4) THEN
	STOP 'Equatorial symmetry is not allowed for the x-boundary'
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Flip signs, y inner boundary !
IF(boundary_flag(3) == 2) THEN
	efac_yin(ixy) = -1
	efac_yin(ixz) = -1
	efac_yin(izx) = -1
	efac_yin(izy) = -1
	ecell_yin(ix) = -1
	ecell_yin(iz) = -1
ELSEIF(boundary_flag(3) == 3) THEN
	efac_yin(iyx) = -1
	efac_yin(iyz) = -1
	efac_yin(izx) = -1
	efac_yin(izy) = -1
	ecell_yin(iy) = -1
	ecell_yin(iz) = -1
ELSEIF(boundary_flag(3) == 4) THEN
	efac_yin(iyx) = -1
	efac_yin(iyz) = -1
	ecell_yin(iy) = -1
END IF

! Flip signs, y outer boundary !
IF(boundary_flag(4) == 2) THEN
	efac_yout(ixy) = -1
	efac_yout(ixz) = -1
	efac_yout(izx) = -1
	efac_yout(izy) = -1
	ecell_yout(ix) = -1
	ecell_yout(iz) = -1
ELSEIF(boundary_flag(4) == 3) THEN
	efac_yout(iyx) = -1
	efac_yout(iyz) = -1
	efac_yout(izx) = -1
	efac_yout(izy) = -1
	ecell_yout(iy) = -1
	ecell_yout(iz) = -1
ELSEIF(boundary_flag(4) == 4) THEN
	efac_yout(iyx) = -1
	efac_yout(iyz) = -1
	ecell_yout(iy) = -1
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Flip signs, z inner boundary !
IF(boundary_flag(5) == 2) THEN
	efac_zin(ixy) = -1
	efac_zin(ixz) = -1
	efac_zin(iyx) = -1
	efac_zin(iyz) = -1
	ecell_zin(ix) = -1
	ecell_zin(iy) = -1
ELSEIF(boundary_flag(5) == 3) THEN
	STOP 'Axial symmetry is not allowed for the z-boundary'
ELSEIF(boundary_flag(5) == 4) THEN
	efac_zin(izx) = -1
	efac_zin(izy) = -1
	ecell_zin(iz) = -1
END IF

! Flip signs, z outer boundary !
IF(boundary_flag(6) == 2) THEN
	efac_zout(ixy) = -1
	efac_zout(ixz) = -1
	efac_zout(iyx) = -1
	efac_zout(iyz) = -1
	ecell_zout(ix) = -1
	ecell_zout(iy) = -1
ELSEIF(boundary_flag(6) == 3) THEN
	STOP 'Axial symmetry is not allowed for the z-boundary'
ELSEIF(boundary_flag(6) == 4) THEN
	efac_zout(izx) = -1
	efac_zout(izy) = -1
	ecell_zout(iz) = -1
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for printing out !

WRITE(*,*)
WRITE(*,*) 'Finished setting up equation indexes'
WRITE(*,*) 'There are', no_of_eq, 'number of equations'
WRITE(*,*)

END SUBROUTINE SETUP_EQN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine allocate arrays related to hydrodynamics variables !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BUILD_HYDRO
USE DEFINITION
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Core arrays for solving hyperbolic PDE

! Allocate arrays for NM !
ALLOCATE(l2(imin2:imax2,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(prim2(imin2:imax2,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(cons2(imin2:imax2,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(u_old2(imin2:imax2,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(prim2_a(imin2:imax2,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))

! Allocate flux arrays for NM !
ALLOCATE(sc2(imin2:imax2,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(flux_2(imin2:imax2,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(dflux_2(imin2:imax2,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For Grid variables !
ALLOCATE (x2(-2:nx_2+3))
ALLOCATE (y2(-2:ny_2+3))
ALLOCATE (z2(-2:nz_2+3))
ALLOCATE (xF2(-2:nx_2+3))
ALLOCATE (yF2(-2:ny_2+3))
ALLOCATE (zF2(-2:nz_2+3))
ALLOCATE (dx2(-2:nx_2+3))
ALLOCATE (dy2(-2:ny_2+3))
ALLOCATE (dz2(-2:nz_2+3))
ALLOCATE (x2bar(-2:nx_2+3))
ALLOCATE (x2cen(-2:nx_2+3))
ALLOCATE (y2cen(-2:ny_2+3))
ALLOCATE (dx2_sq(-2:nx_2+3))
ALLOCATE (dx2_cb(-2:nx_2+3))
ALLOCATE (dx2_qd(-2:nx_2+3))
ALLOCATE (dsin2(-2:ny_2+3))
ALLOCATE (dcos2(-2:ny_2+3))
ALLOCATE (vol2(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))

! Internal energy !
ALLOCATE (epsilon2(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE (eps2_a(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))

! Hydrodynamic variables !
ALLOCATE (cs2(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE (dpdrho2(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE (dpdeps2(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reconstruction weight !

! X-direction !
ALLOCATE (lxm2(-2:nx_2+3))
ALLOCATE (lxm1(-2:nx_2+3))
ALLOCATE (lxc(-2:nx_2+3))
ALLOCATE (lxp1(-2:nx_2+3))
ALLOCATE (lxp2(-2:nx_2+3))
ALLOCATE (rxm2(-2:nx_2+3))
ALLOCATE (rxm1(-2:nx_2+3))
ALLOCATE (rxc(-2:nx_2+3))
ALLOCATE (rxp1(-2:nx_2+3))
ALLOCATE (rxp2(-2:nx_2+3))
ALLOCATE (hpx(-2:nx_2+3))
ALLOCATE (hmx(-2:nx_2+3))

! Y-direction !
ALLOCATE (lym2(-2:ny_2+3))
ALLOCATE (lym1(-2:ny_2+3))
ALLOCATE (lyc(-2:ny_2+3))
ALLOCATE (lyp1(-2:ny_2+3))
ALLOCATE (lyp2(-2:ny_2+3))
ALLOCATE (rym2(-2:ny_2+3))
ALLOCATE (rym1(-2:ny_2+3))
ALLOCATE (ryc(-2:ny_2+3))
ALLOCATE (ryp1(-2:ny_2+3))
ALLOCATE (ryp2(-2:ny_2+3))
ALLOCATE (hpy(-2:ny_2+3))
ALLOCATE (hmy(-2:ny_2+3))

! Z-direction !
ALLOCATE (lzm2(-2:nz_2+3))
ALLOCATE (lzm1(-2:nz_2+3))
ALLOCATE (lzc(-2:nz_2+3))
ALLOCATE (lzp1(-2:nz_2+3))
ALLOCATE (lzp2(-2:nz_2+3))
ALLOCATE (rzm2(-2:nz_2+3))
ALLOCATE (rzm1(-2:nz_2+3))
ALLOCATE (rzc(-2:nz_2+3))
ALLOCATE (rzp1(-2:nz_2+3))
ALLOCATE (rzp2(-2:nz_2+3))
ALLOCATE (hpz(-2:nz_2+3))
ALLOCATE (hmz(-2:nz_2+3))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE