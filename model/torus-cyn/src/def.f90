!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Custom modular files to contain all custom arrays 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE CUSTOM_DEF
USE DEFINITION
IMPLICIT NONE

! floor density, pressure, and epsilon !
REAL*8, ALLOCATABLE, DIMENSION (:,:,:) :: rho_floor
REAL*8, ALLOCATABLE, DIMENSION (:,:,:) :: eps_floor
REAL*8, ALLOCATABLE, DIMENSION (:,:,:) :: p_floor

END MODULE