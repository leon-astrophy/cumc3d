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

! Dummy variables
INTEGER :: i, j, k, l

! Do the inner boundary
IF(boundary_flag(1) == 0) THEN
  DO CONCURRENT(j = 1:3, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2)
    prim2(ivel2_x,nx_min_2-j,k,l) = MIN(prim2(ivel2_x,nx_min_2-j,k,l), 0.0D0)
    prim2(iby,nx_min_2-j,k,l) = 0.0d0
    prim2(ibz,nx_min_2-j,k,l) = 0.0d0             
  ENDDO
ENDIF

! Do the outer boundary
IF(boundary_flag(2) == 0) THEN
  DO CONCURRENT(j = 1:3, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2)
    prim2(ivel2_x,nx_part_2+j,k,l) = MAX(prim2(ivel2_x,nx_part_2+j,k,l), 0.0D0)
    prim2(iby,nx_part_2+j,k,l) = 0.0d0
    prim2(ibz,nx_part_2+j,k,l) = 0.0d0
  ENDDO
ENDIF

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

! Threshold for atmosphere density
REAL*8 :: dphidr
REAL*8 :: rho_min1, rho_min2, factor, diff

! Set up the threshold
rho_min2 = 1.1D0 * prim2_a(irho2)

! Add black hole gravity !
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(factor, diff, dphidr)
DO l = nz_min_2, nz_part_2
  DO k = ny_min_2, ny_part_2
    DO j = nx_min_2, nx_part_2
			diff = prim2(irho2,j,k,l) - rho_min2
      factor = MAX(SIGN(1.0D0, diff), 0.0D0)
      dphidr = 1.0d0/(x2(j) - r_sh)**2
      sc2(ivel2_x,j,k,l) = sc2(ivel2_x,j,k,l) + (-factor*prim2(irho2,j,k,l)*dphidr)
      sc2(itau2,j,k,l) = sc2(itau2,j,k,l) + (-factor*prim2(irho2,j,k,l)*prim2(ivel2_x,j,k,l)*dphidr)
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
