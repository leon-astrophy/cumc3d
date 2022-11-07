!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Written by Leung Shing Chi in 2016
! The subroutine setup the primitive and conservative equation indexes
! Notice that if you want your to add your own quantities, you need to specify them here
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SETUP_EQN
USE DEFINITION
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

! Section for DM !
IF(runDM_flag) THEN

    imin1 = 1

    ! DM density
    no_of_eq = no_of_eq + 1
    imax1 = no_of_eq
    irho1 = no_of_eq
    bfac_x(no_of_eq) = 1
    bfac_y(no_of_eq) = 1
    bfac_z(no_of_eq) = 1
    write(*,*) 'Make irho1 = ', no_of_eq
   
    ! DM vel-x
    no_of_eq = no_of_eq + 1
    imax1 = no_of_eq
    ivel1_x = no_of_eq 
    bfac_x(no_of_eq) = -1
    bfac_y(no_of_eq) = 1
    bfac_z(no_of_eq) = 1
    write(*,*) 'Make ivel1_x = ', no_of_eq

    ! DM vel-y
    no_of_eq = no_of_eq + 1
    imax1 = no_of_eq
    ivel1_y = no_of_eq
    bfac_x(no_of_eq) = 1
    bfac_y(no_of_eq) = -1
    bfac_z(no_of_eq) = 1
    write(*,*) 'Make ivel1_y = ', no_of_eq

    ! DM vel-z
    no_of_eq = no_of_eq + 1
   	imax1 = no_of_eq
   	ivel1_z = no_of_eq
    bfac_x(no_of_eq) = 1
    bfac_y(no_of_eq) = 1
    bfac_z(no_of_eq) = -1
   	write(*,*) 'Make ivel1_z = ', no_of_eq

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now we do the NM sector

! Set up minimum equation for NM !
imin2 = no_of_eq + 1

! NM density
no_of_eq = no_of_eq + 1
imax2 = no_of_eq
irho2 = no_of_eq
bfac_x(no_of_eq) = 1
bfac_y(no_of_eq) = 1
bfac_z(no_of_eq) = 1
write(*,*) 'Make irho2 = ', no_of_eq

! NM vel-x
no_of_eq = no_of_eq + 1
imax2 = no_of_eq
ivel2_x = no_of_eq
bfac_x(no_of_eq) = -1
bfac_y(no_of_eq) = 1
bfac_z(no_of_eq) = 1
write(*,*) 'Make ivel2_x = ', no_of_eq

! NM vel-y
no_of_eq = no_of_eq + 1
imax2 = no_of_eq
ivel2_y = no_of_eq
bfac_x(no_of_eq) = 1
bfac_y(no_of_eq) = -1
bfac_z(no_of_eq) = 1
write(*,*) 'Make ivel2_y = ', no_of_eq

! NM vel-z
no_of_eq = no_of_eq + 1
imax2 = no_of_eq
ivel2_z = no_of_eq
bfac_x(no_of_eq) = 1
bfac_y(no_of_eq) = 1
bfac_z(no_of_eq) = -1
write(*,*) 'Make ivel2_z = ', no_of_eq

! NM epsilon
no_of_eq = no_of_eq + 1
imax2 = no_of_eq
itau2 = no_of_eq
bfac_x(no_of_eq) = 1
bfac_y(no_of_eq) = 1
bfac_z(no_of_eq) = 1
write(*,*) 'Make itau2 = ', no_of_eq

! Dual energies !
IF(dual_energy) THEN
	no_of_eq = no_of_eq + 1
	imax2 = no_of_eq
	ieps2 = no_of_eq
    bfac_x(no_of_eq) = 1
    bfac_y(no_of_eq) = 1
    bfac_z(no_of_eq) = 1
	write(*,*) 'Make ieps2 = ', no_of_eq
END IF

! Dual energies !
IF(etran_flag) THEN
	no_of_eq = no_of_eq + 1
	imax2 = no_of_eq
	iye2 = no_of_eq
    bfac_x(no_of_eq) = 1
    bfac_y(no_of_eq) = 1
    bfac_z(no_of_eq) = 1
	write(*,*) 'Make iye2 = ', no_of_eq
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Core arrays for solving hyperbolic PDE

! Allocate arrays for NM !
ALLOCATE(prim2_a(imin2:imax2))
ALLOCATE(l2(imin2:imax2,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(u2_nm(imin2:imax2,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(u3_nm(imin2:imax2,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(l3_nm(imin2:imax2,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(u_old2(imin2:imax2,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(prim2(imin2:imax2,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(cons2(imin2:imax2,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))

! Allocate flux arrays for NM !
ALLOCATE(sa2(imin2:imax2,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(sb2(imin2:imax2,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(dflux_2(imin2:imax2,1:n_dim,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(flux_2(imin2:imax2,1:n_dim,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))

! For DM !
IF(RUNDM_flag) THEN
    ALLOCATE(prim1_a(imin1:imax1))
    ALLOCATE(l1(min1:imax1,-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
    ALLOCATE(u2_dm(min1:imax1,-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
    ALLOCATE(u3_dm(min1:imax1,-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
    ALLOCATE(l3_dm(min1:imax1,-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
    ALLOCATE(u_old1(min1:imax1,-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
    ALLOCATE(prim1(min1:imax1,-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
    ALLOCATE(cons1(min1:imax1,-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
    ALLOCATE(sa1(min1:imax1,-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
	ALLOCATE(sb1(min1:imax1,-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
	ALLOCATE(flux_1(imin1:imax1,1:n_dim,2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
	ALLOCATE(dflux_1(imin1:imax1,1:n_dim,2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
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

! Legendre polynominal and mutlipole moments !
IF(coordinate_flag == 2) THEN
	ALLOCATE (inner_2(0:2*lmax,-2:nx_2+3))
	ALLOCATE (outer_2(0:2*lmax,-2:nx_2+3))
	ALLOCATE (rsolid_2(0:2*lmax,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
	ALLOCATE (isolid_2(0:2*lmax,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ELSE
    ALLOCATE (qpole2(0:2*lmax))
END IF
ALLOCATE (legendre2(0:2*lmax,2:nx_2+3,-2:ny_2+3,-2:nz_2+3))

! For Grid variables !
ALLOCATE (dr2(-2:nx_2+3))
ALLOCATE (x2(-2:nx_2+3))
ALLOCATE (y2(-2:ny_2+3))
ALLOCATE (z2(-2:nz_2+3))
ALLOCATE (xF2(-2:nx_2+3))
ALLOCATE (yF2(-2:ny_2+3))
ALLOCATE (zF2(-2:nz_2+3))
ALLOCATE (rad2(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE (cos2(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE (sin2(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE (vol2(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE (radbar2(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE (volbar2(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))

! Internal energy !
ALLOCATE (epsilon2(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))

! Hydrodynamic variables !
ALLOCATE (cs2(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE (dpdrho2(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE (dpdeps2(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))

! Temperature !
IF(have_temp) THEN
	ALLOCATE (temp2(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
END IF

! Potentials !
ALLOCATE (phi2(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE (phi2_dm(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE (phi2_nm(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE (dphi_2(1:n_dim,2:nx_2+3,-2:ny_2+3,-2:nz_2+3))

! Extra Hydrodynamic variables !
IF (dual_energy) THEN
    ALLOCATE (dp_2(1:n_dim,2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For the DM sector !

IF (DM_flag) THEN

	! Legendre polynominal and mutlipole moments !
	IF(coordinate_flag == 2) THEN
		ALLOCATE (rsolid_1(0:2*lmax,2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
		ALLOCATE (isolid_1(0:2*lmax,2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
		ALLOCATE (inner_1(0:2*lmax,-2:nx_1+3))
		ALLOCATE (outer_1(0:2*lmax,-2:nx_1+3))
	ELSE
        ALLOCATE (qpole1(0:2*lmax))
    END IF
    ALLOCATE (legendre1(0:2*lmax,-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))

	! For Grid variables !
	ALLOCATE (dr1(-2:nx_1+3))
	ALLOCATE (x1(-2:nx_1+3))
	ALLOCATE (y1(-2:ny_1+3))
	ALLOCATE (z1(-2:nz_1+3))
	ALLOCATE (xF1(-2:nx_1+3))
	ALLOCATE (yF1(-2:ny_1+3))
	ALLOCATE (zF1(-2:nz_1++3))
	ALLOCATE (rad1(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
	ALLOCATE (cos1(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
	ALLOCATE (sin1(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
	ALLOCATE (vol1(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
	ALLOCATE (radbar1(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
	ALLOCATE (volbar1(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))

	! Hydrodynamic variables !
	ALLOCATE (p1(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
	ALLOCATE (cs1(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
	ALLOCATE (dpdrho1(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
	ALLOCATE (dpdeps1(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))

	! Potentials !
	ALLOCATE (phi1(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
	ALLOCATE (phi1_dm(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
	ALLOCATE (phi1_nm(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
	ALLOCATE (dphi_1(1:n_dim,2:nx_1+3,-2:ny_1+3,-2:nz_1+3))

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE
