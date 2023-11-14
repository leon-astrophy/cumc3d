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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Equations identifiers !

! Minimum/Maximum Eq. number for NM
INTEGER :: imin
INTEGER :: imax

! Identifiers for the NM variables !
! NM density !
INTEGER :: irho

! NM x-velocity
INTEGER :: ivx

! NM y-velocity
INTEGER :: ivy !so 

! NM z-velocity
INTEGER :: ivz

! NM total energy density
INTEGER :: itau

! magnetic fields !
INTEGER :: ibx, iby, ibz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\
! Grid variables !

! R and Z coordinate of the grid for NM !
REAL*8, ALLOCATABLE, DIMENSION (:) :: x
REAL*8, ALLOCATABLE, DIMENSION (:) :: y
REAL*8, ALLOCATABLE, DIMENSION (:) :: z
REAL*8, ALLOCATABLE, DIMENSION (:) :: xF
REAL*8, ALLOCATABLE, DIMENSION (:) :: yF
REAL*8, ALLOCATABLE, DIMENSION (:) :: zF
REAL*8, ALLOCATABLE, DIMENSION (:) :: dx
REAL*8, ALLOCATABLE, DIMENSION (:) :: dy
REAL*8, ALLOCATABLE, DIMENSION (:) :: dz
REAL*8, ALLOCATABLE, DIMENSION (:) :: xbar
REAL*8, ALLOCATABLE, DIMENSION (:) :: dx_cb
REAL*8, ALLOCATABLE, DIMENSION (:) :: sine
REAL*8, ALLOCATABLE, DIMENSION (:) :: sinf
REAL*8, ALLOCATABLE, DIMENSION (:) :: dsine
REAL*8, ALLOCATABLE, DIMENSION (:) :: dcose
REAL*8, ALLOCATABLE, DIMENSION (:,:,:) :: vol

! Reconstruction weight !
REAL*8, ALLOCATABLE, DIMENSION (:,:) :: wx
REAL*8, ALLOCATABLE, DIMENSION (:,:) :: wy
REAL*8, ALLOCATABLE, DIMENSION (:,:) :: wz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The following are the hydro set for NM sector

! NM primitive and conservative variables !
REAL*8, ALLOCATABLE, DIMENSION (:,:,:,:) :: prim
REAL*8, ALLOCATABLE, DIMENSION (:,:,:,:) :: cons

! NM speed of sound, and pressure derivatives !
REAL*8, ALLOCATABLE, DIMENSION (:,:,:) :: cs
REAL*8, ALLOCATABLE, DIMENSION (:,:,:) :: dpdrho
REAL*8, ALLOCATABLE, DIMENSION (:,:,:) :: dpdeps

! NM internal energy !
REAL*8, ALLOCATABLE, DIMENSION (:,:,:) :: epsilon

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for time-evolution !

! For RK-Time evolution 
REAL*8, ALLOCATABLE, DIMENSION (:,:,:,:) :: l_rk, u_old

! Flux arrays for NM !
REAL*8, allocatable, DIMENSION (:,:,:,:) :: sc
REAL*8, allocatable, DIMENSION (:,:,:,:) :: flux

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Real/integer variable scalars !

! Time step !
REAL*8 :: dt

! RK step number !
INTEGER :: n_step

! Iteration step number
INTEGER :: n_iter

! Polytropic index !
REAL*8 :: kgas, ggas

! for openacc !
#ifdef GPU
!$acc declare create (kgas)
!$acc declare create (ggas)
#endif

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