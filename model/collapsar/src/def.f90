!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Custom modular files to contain all custom arrays 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE CUSTOM_DEF
USE DEFINITION
IMPLICIT NONE

! gravitational potential !
REAL*8, ALLOCATABLE, DIMENSION (:,:,:) :: phi, phi_old

! for poisson solver !
REAL*8, ALLOCATABLE, DIMENSION (:) :: ajp1, ajm1
REAL*8, ALLOCATABLE, DIMENSION (:,:) :: bkp1, bkm1
REAL*8, ALLOCATABLE, DIMENSION (:,:,:) :: clp1, clm1, epsc

! for poisson solver !
REAL*8, PARAMETER :: omega_weight = 1.9d0

! mass of black hole !
REAL*8 :: m_bh, mbh_old

! mass injection due to floor !
REAL*8 :: m_inj

END MODULE
