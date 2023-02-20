!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! This module contains all the arrays and variables that are neccessary to run the hydro
! simulations of either pure hydro, discs, or self gravitating fluid objects like stars
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE DEFINITION
IMPLICIT NONE
SAVE
INCLUDE "parameter.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Boundary flag for variables !

INTEGER :: bfac_xin(50)
INTEGER :: bfac_yin(50)
INTEGER :: bfac_zin(50)
INTEGER :: bfac_xout(50)
INTEGER :: bfac_yout(50)
INTEGER :: bfac_zout(50)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Equations identifiers !

! Minimum/Maximum Eq. number for DM
INTEGER :: imin1
INTEGER :: imax1

! Identifiers for the DM variables !

! DM density !
INTEGER :: irho1

! DM x-velocity
INTEGER :: ivel1_x

! DM y-velocity
INTEGER :: ivel1_y

! DM z-velocity
INTEGER :: ivel1_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Equations identifiers !

! Minimum/Maximum Eq. number for NM
INTEGER :: imin2
INTEGER :: imax2

! Identifiers for the NM variables !

! NM density !
INTEGER :: irho2

! NM x-velocity
INTEGER :: ivel2_x

! NM y-velocity
INTEGER :: ivel2_y

! NM z-velocity
INTEGER :: ivel2_z

! NM total energy density
INTEGER :: itau2

! NM internal energy density
INTEGER :: ieps2

! NM electron fractions
INTEGER :: iye2

! magnetic fields !
INTEGER :: ibx, iby, ibz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\
! Grid variables !

! Grid coordinates for DM !
REAL*8, ALLOCATABLE, DIMENSION (:) :: x1
REAL*8, ALLOCATABLE, DIMENSION (:) :: y1
REAL*8, ALLOCATABLE, DIMENSION (:) :: z1
REAL*8, ALLOCATABLE, DIMENSION (:) :: xF1
REAL*8, ALLOCATABLE, DIMENSION (:) :: yF1
REAL*8, ALLOCATABLE, DIMENSION (:) :: zF1
REAL*8, ALLOCATABLE, DIMENSION (:) :: dx1
REAL*8, ALLOCATABLE, DIMENSION (:) :: dy1
REAL*8, ALLOCATABLE, DIMENSION (:) :: dz1
REAL*8, ALLOCATABLE, DIMENSION (:,:,:) :: vol1

! R and Z coordinate of the grid for NM !
REAL*8, ALLOCATABLE, DIMENSION (:) :: x2
REAL*8, ALLOCATABLE, DIMENSION (:) :: y2
REAL*8, ALLOCATABLE, DIMENSION (:) :: z2
REAL*8, ALLOCATABLE, DIMENSION (:) :: xF2
REAL*8, ALLOCATABLE, DIMENSION (:) :: yF2
REAL*8, ALLOCATABLE, DIMENSION (:) :: zF2
REAL*8, ALLOCATABLE, DIMENSION (:) :: dx2
REAL*8, ALLOCATABLE, DIMENSION (:) :: dy2
REAL*8, ALLOCATABLE, DIMENSION (:) :: dz2
REAL*8, ALLOCATABLE, DIMENSION (:,:,:) :: vol2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Hydrodynamical variables

! Geometric factor for cylindrical/spherical coordiantes !
REAL*8, ALLOCATABLE, DIMENSION (:) :: geom1_x, geom1_y

! DM atmospheric primitive variables !
REAL*8, ALLOCATABLE, DIMENSION (:) :: prim1_a

! DM primitive and conservative variables !
REAL*8, ALLOCATABLE, DIMENSION (:,:,:,:) :: prim1
REAL*8, ALLOCATABLE, DIMENSION (:,:,:,:) :: cons1

! DM pressure, speed of sound, and pressure derivatives !
REAL*8, ALLOCATABLE, DIMENSION (:,:,:) :: p1
REAL*8, ALLOCATABLE, DIMENSION (:,:,:) :: cs1
REAL*8, ALLOCATABLE, DIMENSION (:,:,:) :: dpdrho1
REAL*8, ALLOCATABLE, DIMENSION (:,:,:) :: dpdeps1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The following are the hydro set for NM sector

! Geometric factor for cylindrical/spherical coordiantes !
REAL*8, ALLOCATABLE, DIMENSION (:) :: geom2_x, geom2_y

! NM atmospheric primitive variables !
REAL*8, ALLOCATABLE, DIMENSION (:) :: prim2_a

! NM primitive and conservative variables !
REAL*8, ALLOCATABLE, DIMENSION (:,:,:,:) :: prim2
REAL*8, ALLOCATABLE, DIMENSION (:,:,:,:) :: cons2

! NM speed of sound, and pressure derivatives !
REAL*8, ALLOCATABLE, DIMENSION (:,:,:) :: cs2
REAL*8, ALLOCATABLE, DIMENSION (:,:,:) :: dpdrho2
REAL*8, ALLOCATABLE, DIMENSION (:,:,:) :: dpdeps2

! NM internal energy !
REAL*8, ALLOCATABLE, DIMENSION (:,:,:) :: epsilon2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for time-evolution !

! For RK-Time evolution 
REAL*8, ALLOCATABLE, DIMENSION (:,:,:,:) :: l1, u_old1
REAL*8, ALLOCATABLE, DIMENSION (:,:,:,:) :: l2, u_old2

! The auxillary array for the flux term, DM
REAL*8, allocatable, DIMENSION (:,:,:,:) :: sc1
REAL*8, allocatable, DIMENSION (:,:,:,:) :: flux_1
REAL*8, allocatable, DIMENSION (:,:,:,:) :: dflux_1

! Flux arrays for NM !
REAL*8, allocatable, DIMENSION (:,:,:,:) :: sc2
REAL*8, allocatable, DIMENSION (:,:,:,:) :: flux_2
REAL*8, allocatable, DIMENSION (:,:,:,:) :: dflux_2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Real/integer variable scalars !

! The number of highest grid along three axis !
INTEGER :: nx_part_1, ny_part_1, nz_part_1
INTEGER :: nx_part_2, ny_part_2, nz_part_2
							
! The number of lowest grid along three axis !
INTEGER :: nx_min_1, ny_min_1, nz_min_1
INTEGER :: nx_min_2, ny_min_2, nz_min_2

! Time step !
REAL*8 :: dt

! Iteration step number
INTEGER :: n_iter

! Polytropic index !
REAL*8 :: kgas1, ggas1
REAL*8 :: kgas2, ggas2

! Atmospheric epsilon !
REAL*8 :: eps1_a, eps2_a

! Global simulation time !
REAL*8 :: global_time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Storing coefficients for the RK method !
REAL*8, PARAMETER :: rk20 = 3.0d0/4.0d0
REAL*8, PARAMETER :: rk21 = 1.0d0/4.0d0
REAL*8, PARAMETER :: rk22 = 1.0d0/4.0d0
REAL*8, PARAMETER :: rk30 = 1.0d0/3.0d0
REAL*8, PARAMETER :: rk31 = 2.0d0/3.0d0
REAL*8, PARAMETER :: rk32 = 2.0d0/3.0d0

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! function for computing cartesian dot_product
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	REAL*8 function dot_product(vec_a, vec_b)
	implicit none
	REAL*8, DIMENSION (1:3) :: vec_a, vec_b
	dot_product = vec_a(1)*vec_b(1) + vec_a(2)*vec_b(2) + vec_a(3)*vec_b(3)
	end function

END MODULE DEFINITION