!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom equations !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_EQN
USE DEFINITION
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom arrays !
!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_HYDRO
USE CUSTOM_DEF
USE DEFINITION
IMPLICIT NONE

! gravitational potential energy !
ALLOCATE (rho_floor(-2:nx+3,-2:ny+3,-2:nz+3))
ALLOCATE (p_floor(-2:nx+3,-2:ny+3,-2:nz+3))
ALLOCATE (eps_floor(-2:nx+3,-2:ny+3,-2:nz+3))

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Populate custom arrays to GPU !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_POPULATE
USE DEFINITION
USE CUSTOM_DEF
IMPLICIT NONE

! Now populate all necessary, and reuseable arrays to the graphic cards !
!$ACC enter DATA COPYIN(rho_floor, p_floor, eps_floor)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Clear custom arrays from GPU !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_CLEAR
USE DEFINITION
USE CUSTOM_DEF
IMPLICIT NONE

! Now we clear memory in the GPU device !
!$ACC exit DATA DELETE(rho_floor, p_floor, eps_floor)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom grid !
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_GRID
USE DEFINITION
IMPLICIT NONE
INCLUDE "param.h"

! Integer !
INTEGER :: i, nlines
INTEGER :: j, k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(refine_grid) THEN
  ! Read the number of lines in the file !
  nlines = 0 
  OPEN (999, file = './grid/r_grid.dat') 
  DO 
    READ (999,*, END=10) 
    nlines = nlines + 1 
  END DO 
  10 CLOSE (999) 

  ! Error message !
  IF(nlines .ne. nx+7) THEN
    WRITE (*,*) 'number of r-grid faces from files', nlines-7
    WRITE (*,*) 'number of r-grid faces in the program', nx
    STOP 'inconsistent number of r-grid faces, exit'
  END IF

  ! Read !
  OPEN(UNIT=999, FILE = './grid/r_grid.dat', ACTION='READ')
  DO i = -3, nx+3
    READ(999,*) xF(i)
  ENDDO
  CLOSE(999)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Read the number of lines in the file !
  nlines = 0 
  OPEN (999, file = './grid/phi_grid.dat') 
  DO 
    READ (999,*, END=20) 
    nlines = nlines + 1 
  END DO 
  20 CLOSE (999) 

  ! Error message !
  IF(nlines .ne. ny+7) THEN
    WRITE (*,*) 'number of phi-grid faces from files', nlines-7
    WRITE (*,*) 'number of phi-grid faces in the program', ny
    STOP 'inconsistent number of phi-grid faces, exit'
  END IF

  ! Read !
  OPEN(UNIT=999, FILE = './grid/phi_grid.dat', ACTION='READ')
  DO i = -3, ny+3
    READ(999,*) yF(i)
  ENDDO
  CLOSE(999)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Read the number of lines in the file !
  nlines = 0 
  OPEN (999, file = './grid/z_grid.dat') 
  DO 
    READ (999,*, END=30) 
    nlines = nlines + 1 
  END DO 
  30 CLOSE (999) 

  ! Error message !
  IF(nlines .ne. nz+7) THEN
    WRITE (*,*) 'number of z-grid faces from files', nlines-7
    WRITE (*,*) 'number of z-grid faces in the program', nz
    STOP 'inconsistent number of z-grid faces, exit'
  END IF

  ! Read !
  OPEN(UNIT=999, FILE = './grid/z_grid.dat', ACTION='READ')
  DO i = -3, nz+3
    READ(999,*) zF(i)
  ENDDO
  CLOSE(999)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom grid !
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_BOUNDARY_X
USE DEFINITION
IMPLICIT NONE

! Dummy variables
INTEGER :: i, j, k, l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT)   
DO l = -2, nz + 3
  DO k = -2, ny + 3
    DO j = 1, 3
      prim(ivx,1-j,k,l) = MIN(prim(ivx,1-j,k,l), 0.0D0)
      prim(ivx,nx+j,k,l) = MAX(prim(ivx,nx+j,k,l), 0.0D0)
    END DO
  END DO               
ENDDO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom grid !
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_BOUNDARY_Y
USE DEFINITION
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom grid !
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_BOUNDARY_Z
USE DEFINITION
IMPLICIT NONE

! Dummy variables
INTEGER :: i, j, k, l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT)  
DO l = 1, 3
  DO k = -2, ny + 3
    DO j = -2, nx + 3
      prim(ivz,j,k,1-l) = MIN(prim(ivz,j,k,1-l), 0.0d0)
      prim(ivz,j,k,nz+l) = MAX(prim(ivz,j,k,nz+l), 0.0d0)
    END DO
  END DO               
ENDDO 
!$ACC END PARALLEL
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Back up fluxes from riemann solvers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETFLUX_X
USE DEFINITION
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Back up fluxes from riemann solvers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETFLUX_Y
USE DEFINITION
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Back up fluxes from riemann solvers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETFLUX_Z
USE DEFINITION
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!
! Custom variable floor !
!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_CHECKRHO
USE CUSTOM_DEF
USE DEFINITION
IMPLICIT NONE
  
! Dummy variables
INTEGER :: i, j, k, l

! Threshold for atmosphere density
REAL*8 :: factor, diff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP PARALLEL DO COLLAPSE(3) PRIVATE(diff, factor) SCHEDULE(STATIC) 
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) PRIVATE(diff, factor)
DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx

      ! Check density !
      diff = prim(irho,j,k,l) - rho_floor(j,k,l)
      factor = MAX(SIGN(1.0D0, diff), 0.0D0)
      prim(irho,j,k,l) = factor*prim(irho,j,k,l) + (1.0D0 - factor)*rho_floor(j,k,l)
      prim(ivx:ivz,j,k,l) = factor*prim(ivx:ivz,j,k,l)
      epsilon(j,k,l) = factor*epsilon(j,k,l) + (1.0D0 - factor)*eps_floor(j,k,l)
      IF(epsilon(j,k,l) < 0.0d0) epsilon(j,k,l) = eps_floor(j,k,l)
      
    END DO
  END DO
END DO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom equations !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_SOURCE
USE DEFINITION
USE CUSTOM_DEF
IMPLICIT NONE
INCLUDE "param.h"

! Integer !
INTEGER :: i, j, k, l

! real !
REAL*8 :: radius

! Threshold for atmosphere density
REAL*8 :: dphidr, dphidz
REAL*8 :: factor, diff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Add black hole gravity !
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(radius, dphidr, dphidz, factor, diff)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) PRIVATE(radius, dphidr, dphidz, factor, diff)
DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx
      diff = prim(irho,j,k,l) - rho_floor(j,k,l)
      factor = MAX(SIGN(1.0D0, diff), 0.0D0)
      radius = DSQRT(x(j)*x(j) + z(l)*z(l))
      dphidr = x(j)/radius/((radius - r_sh)*(radius - r_sh))
      dphidz = z(l)/radius/((radius - r_sh)*(radius - r_sh))
      sc(ivx,j,k,l) = sc(ivx,j,k,l) + (-factor*prim(irho,j,k,l)*dphidr)
      sc(ivz,j,k,l) = sc(ivz,j,k,l) + (-factor*prim(irho,j,k,l)*dphidz)
      sc(itau,j,k,l) = sc(itau,j,k,l) + (-factor*prim(irho,j,k,l)*(prim(ivx,j,k,l)*dphidr + prim(ivz,j,k,l)*dphidz))
    END DO
  END DO
END DO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!
! Do custom updates !
!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_UPDATE (p_in)
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER, INTENT(IN) :: p_in

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom equations !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OPERATOR_SPLIT
USE DEFINITION
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom grid !
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OPENFILE_CUSTOM
USE DEFINITION
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom grid !
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_ANALYSIS
USE DEFINITION
IMPLICIT NONE

END SUBROUTINE
