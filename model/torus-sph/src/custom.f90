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
USE DEFINITION
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Populate custom arrays to GPU !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_POPULATE
USE DEFINITION
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Clear custom arrays from GPU !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_CLEAR
USE DEFINITION
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom grid !
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_GRID
USE DEFINITION
IMPLICIT NONE
INCLUDE "param.h"

! Intger !
INTEGER :: i, j, k, l
INTEGER :: j1, j2

! Real !
REAL*8 :: s_grid
REAL*8 :: dx_fine, dy_fine
REAL*8 :: alpha, beta
REAL*8 :: r1, r2

! Do refined grid if neccessary !
IF(refined_grid) THEN

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! x-direction 

  ! Parameter defining the number of grid !
  j1 = x_grid
  j2 = nx_2

  ! Define grid size !
  dx_fine = (x_fine - x2_start)/DBLE(x_grid)

  ! Assign finest resolution grid !
  DO j = -2, x_grid
    xF2(j) = x2_start + DBLE(j)*dx_fine
  END DO

  ! Exponential parameters
  r1 = x_fine
  r2 = x2_end
  beta = 1.0d0/DBLE(j2 - j1)*log(r2/r1)
  alpha = r1*exp(-DBLE(j1)/DBLE(j2-j1)*log(r2/r1))

  ! Exponential grid for the outer domain !
  DO j = x_grid + 1, nx_2 + 3
    xF2(j) = alpha*exp(beta*DBLE(j))
  END DO

END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom grid !
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_BOUNDARY_X
USE DEFINITION
IMPLICIT NONE

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

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!
! Custom variable floor !
!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOMFLOOR
USE DEFINITION
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom equations !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_SOURCE
USE DEFINITION
IMPLICIT NONE
INCLUDE "param.h"

! Integer !
INTEGER :: i, j, k, l

! Threshold for atmosphere density
REAL*8 :: dphidr
REAL*8 :: rho_min1, rho_min2, factor, diff

! Add black hole gravity !
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(factor, diff, dphidr, rho_min2)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) PRIVATE(factor, diff, dphidr, rho_min2)
DO l = nz_min_2, nz_part_2
  DO k = ny_min_2, ny_part_2
    DO j = nx_min_2, nx_part_2
      rho_min2 = 1.1D0 * prim2_a(irho2)
			diff = prim2(irho2,j,k,l) - rho_min2
      factor = MAX(SIGN(1.0D0, diff), 0.0D0)
      dphidr = 1.0d0/(x2(j) - r_sh)**2
      sc2(ivel2_x,j,k,l) = sc2(ivel2_x,j,k,l) + (-factor*prim2(irho2,j,k,l)*dphidr)
      sc2(itau2,j,k,l) = sc2(itau2,j,k,l) + (-factor*prim2(irho2,j,k,l)*prim2(ivel2_x,j,k,l)*dphidr)
    END DO
  END DO
END DO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

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
