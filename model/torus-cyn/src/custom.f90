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

! Set up the threshold
rho_min2 = 1.1D0 * prim2_a(irho2)

! Add black hole gravity !
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(radius, factor, diff, dphidr, dphidz)
DO l = nz_min_2, nz_part_2
  DO k = ny_min_2, ny_part_2
    DO j = nx_min_2, nx_part_2
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
