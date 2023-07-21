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
REAL*8, ALLOCATABLE, DIMENSION (:) :: x1bar
REAL*8, ALLOCATABLE, DIMENSION (:) :: x1cen
REAL*8, ALLOCATABLE, DIMENSION (:) :: y1cen
REAL*8, ALLOCATABLE, DIMENSION (:) :: dx1_sq
REAL*8, ALLOCATABLE, DIMENSION (:) :: dx1_cb
REAL*8, ALLOCATABLE, DIMENSION (:) :: dx1_qd
REAL*8, ALLOCATABLE, DIMENSION (:) :: dsin1
REAL*8, ALLOCATABLE, DIMENSION (:) :: dcos1
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
REAL*8, ALLOCATABLE, DIMENSION (:) :: x2cen
REAL*8, ALLOCATABLE, DIMENSION (:) :: y2cen
REAL*8, ALLOCATABLE, DIMENSION (:) :: x2bar
REAL*8, ALLOCATABLE, DIMENSION (:) :: dx2_sq
REAL*8, ALLOCATABLE, DIMENSION (:) :: dx2_cb
REAL*8, ALLOCATABLE, DIMENSION (:) :: dx2_qd
REAL*8, ALLOCATABLE, DIMENSION (:) :: sin2
REAL*8, ALLOCATABLE, DIMENSION (:) :: sin2f
REAL*8, ALLOCATABLE, DIMENSION (:) :: dsin2
REAL*8, ALLOCATABLE, DIMENSION (:) :: dcos2
REAL*8, ALLOCATABLE, DIMENSION (:,:,:) :: vol2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Hydrodynamical variables

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

! DM internal energy !
REAL*8, ALLOCATABLE, DIMENSION (:,:,:) :: epsilon1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The following are the hydro set for NM sector

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
! Reconstruction weight !

! h-parameter for PPM !
REAL*8, ALLOCATABLE, DIMENSION (:) :: hpx, hpy, hpz
REAL*8, ALLOCATABLE, DIMENSION (:) :: hmx, hmy, hmz

! NM along the x, y, z direction !
REAL*8, ALLOCATABLE, DIMENSION (:) :: lxm2, lxm1, lxc, lxp1, lxp2
REAL*8, ALLOCATABLE, DIMENSION (:) :: lym2, lym1, lyc, lyp1, lyp2
REAL*8, ALLOCATABLE, DIMENSION (:) :: lzm2, lzm1, lzc, lzp1, lzp2
REAL*8, ALLOCATABLE, DIMENSION (:) :: rxm2, rxm1, rxc, rxp1, rxp2
REAL*8, ALLOCATABLE, DIMENSION (:) :: rym2, rym1, ryc, ryp1, ryp2
REAL*8, ALLOCATABLE, DIMENSION (:) :: rzm2, rzm1, rzc, rzp1, rzp2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for time-evolution !

! For RK-Time evolution 
REAL*8, ALLOCATABLE, DIMENSION (:,:,:,:) :: l1, u_old1
REAL*8, ALLOCATABLE, DIMENSION (:,:,:,:) :: l2, u_old2

! The auxillary array for the flux term, DM
REAL*8, allocatable, DIMENSION (:,:,:,:) :: sc1
REAL*8, allocatable, DIMENSION (:,:,:,:) :: flux_1

! Flux arrays for NM !
REAL*8, allocatable, DIMENSION (:,:,:,:) :: sc2
REAL*8, allocatable, DIMENSION (:,:,:,:) :: flux_2

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

! RK step number !
INTEGER :: n_step

! Iteration step number
INTEGER :: n_iter

! Polytropic index !
REAL*8 :: kgas1, ggas1
REAL*8 :: kgas2, ggas2

! for openacc !
#ifdef GPU
!$acc declare create (kgas2)
!$acc declare create (ggas2)
#endif

! Global simulation time !
REAL*8 :: global_time

! Atmospheric epsilon !
REAL*8 :: eps1_a, eps2_a

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

	! dot product between two vectors !
	REAL*8 function dot_product(vec_a, vec_b)
	!$ACC routine seq
	implicit none
	REAL*8, DIMENSION (1:3) :: vec_a, vec_b
	dot_product = vec_a(1)*vec_b(1) + vec_a(2)*vec_b(2) + vec_a(3)*vec_b(3)
	end function

	! the kroncker delta function !
	REAL*8 function kroncker_delta(i, j)
	!$ACC routine seq
	implicit none
	INTEGER :: i, j
		IF(i == j) THEN
			kroncker_delta = 1.0d0
		ELSE
			kroncker_delta = 0.0D0
		END IF
	end function

END MODULE DEFINITION