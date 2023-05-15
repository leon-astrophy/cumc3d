!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom equations !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_EQN
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
REAL*8 :: dx_fine, dz_fine
REAL*8 :: alpha, beta
REAL*8 :: r1, r2

! Do refined grid if neccessary !
IF(refined_grid) THEN

  ! Define grid size !
  dx_fine = x_fine/DBLE(x_grid)

  ! exponential parameters
  r1 = x2_start + x_fine
  r2 = x2_end
  j1 = x_grid
  j2 = nx_2
  beta = 1.0d0/DBLE(j2 - j1)*log(r2/r1)
  alpha = r1*exp(-DBLE(j1)/DBLE(j2-j1)*log(r2/r1))

  ! My custom grid !
  DO j = -2, x_grid
    xF2(j) = x2_start + DBLE(j)*dx_fine
  END DO
  DO j = x_grid + 1, nx_2 + 3
    xF2(j) = alpha*exp(beta*DBLE(j))
  END DO

  ! Define grid size !
  dz_fine = z_fine/DBLE(z_grid)

  ! exponential parmaeters
  r1 = z_fine/2
  r2 = z2_end
  j1 = nz_2/2 + z_grid/2
  j2 = nz_2
  beta = 1.0d0/DBLE(j2 - j1)*log(r2/r1)
  alpha = r1*exp(-DBLE(j1)/DBLE(j2-j1)*log(r2/r1))

  ! My custom grid !
  DO j = nz_2/2, nz_2/2 + z_grid/2
    zF2(j) = DBLE(j - nz_2/2)*dz_fine
  END DO
  DO j = nz_2/2 + z_grid/2, nz_2 + 3
    zF2(j) = alpha*exp(beta*DBLE(j))
  END DO

  ! exponential parmaeters
  r1 = - z_fine/2
  r2 = z2_start
  j1 = nz_2/2 - z_grid/2
  j2 = 1
  beta = 1.0d0/DBLE(j2 - j1)*log(r2/r1)
  alpha = r1*exp(-DBLE(j1)/DBLE(j2-j1)*log(r2/r1))

  ! My custom grid !
  DO j = nz_2/2, nz_2/2 - z_grid/2, -1
    zF2(j) = DBLE(j - nz_2/2)*dz_fine
  END DO
  DO j = nz_2/2 - z_grid/2, -2, -1
    zF2(j) = alpha*exp(beta*DBLE(j))
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

! real !
REAL*8 :: radius

! Threshold for atmosphere density
REAL*8 :: dphidr, dphidz
REAL*8 :: rho_min1, rho_min2, factor, diff

! Add black hole gravity !
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(radius, factor, diff, dphidr, dphidz, rho_min2)
DO l = nz_min_2, nz_part_2
  DO k = ny_min_2, ny_part_2
    DO j = nx_min_2, nx_part_2
      rho_min2 = 1.1D0 * prim2_a(irho2,j,k,l)
      diff = prim2(irho2,j,k,l) - rho_min2
      factor = MAX(SIGN(1.0D0, diff), 0.0D0)
      radius = SQRT(x2(j)**2 + z2(l)**2)
      dphidr = x2(j)/radius/(radius - r_sh)**2
      dphidz = z2(l)/radius/(radius - r_sh)**2
      sc2(ivel2_x,j,k,l) = sc2(ivel2_x,j,k,l) + (-factor*prim2(irho2,j,k,l)*dphidr)
      sc2(ivel2_z,j,k,l) = sc2(ivel2_z,j,k,l) + (-factor*prim2(irho2,j,k,l)*dphidz)
      sc2(itau2,j,k,l) = sc2(itau2,j,k,l) + (-factor*prim2(irho2,j,k,l)*(prim2(ivel2_x,j,k,l)*dphidr + prim2(ivel2_z,j,k,l)*dphidz))
    END DO
  END DO
END DO
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