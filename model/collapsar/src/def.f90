!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Custom modular files to contain all custom arrays 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE CUSTOM_DEF
USE DEFINITION
IMPLICIT NONE

! Unit conversion !
REAL*8 :: vel2code, gauss2code
REAL*8 :: masscgs2code, lengthcgs2code, tcgs2code

! Floors !
REAL*8, ALLOCATABLE, DIMENSION (:) :: prim_a

! Mass flux !
REAL*8, ALLOCATABLE, DIMENSION (:,:) :: mflux_x

! gravitational potential !
REAL*8, ALLOCATABLE, DIMENSION (:,:,:) :: phi, phi_old

! for poisson solver !
REAL*8, ALLOCATABLE, DIMENSION (:) :: ajp1, ajm1
REAL*8, ALLOCATABLE, DIMENSION (:,:) :: bkp1, bkm1
REAL*8, ALLOCATABLE, DIMENSION (:,:,:) :: clp1, clm1, epsc

! for poisson solver !
REAL*8, PARAMETER :: omega_weight = 1.90d0 !6588019511d0

! schwarzschild radius !
REAL*8 :: r_sh

! mass of black hole !
REAL*8 :: m_bh
REAL*8 :: mbh_old

! mass injection due to floor !
REAL*8 :: m_inj

! magnetic flux at horizon !
REAL*8 :: bflux

! epsilon at atmosphere !
REAL*8 :: eps_a

END MODULE
