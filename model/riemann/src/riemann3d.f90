!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! 
! 3D hydro code test
! see. Mon. Not. R. Astron. Soc. 390, 1267–1281 (2008)
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Riemann_3d
USE DEFINITION
IMPLICIT NONE
INCLUDE "param.h"

! Integer and real numbers !
INTEGER :: i, j, k, l
REAL*8 :: dummy

! 3D Explosion with gamma = 1.4
IF(test_model == 1) THEN

  ggas2 = 1.4D0

  DO j = 1, nx_2
    IF(x2(j) < 0.4D0) THEN
      prim2(irho2,j,:,:) = 1.0D0       
      prim2(itau2,j,:,:) = 1.0D0
      prim2(ivel2_x,j,:,:) = 0.0D0
      prim2(ivel2_y,j,:,:) = 0.0D0
      epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2 - 1.0D0)
    ELSE
      prim2(irho2,j,:,:) = 0.125D0
      prim2(itau2,j,:,:) = 0.1D0
      prim2(ivel2_x,j,:,:) = 0.0D0
      prim2(ivel2_y,j,:,:) = 0.0D0
      epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2 - 1.0D0)
    END IF
  END DO

! Oblique Sod shock tube problem
ELSEIF(test_model == 2) THEN

  ggas2 = (5.0D0/3.0D0)

  DO l = 1, nz_2
    DO k = 1, ny_2
      DO j = 1, nx_2
        IF(x2(j) + y2(k) + z2(l) < 0.5D0) THEN
          prim2(irho2,j,k,l) = 4.0D0    
          prim2(itau2,j,k,l) = 1.0D0
        ELSEIF(x2(j) + y2(k) + z2(l) > 0.5D0 .AND. x2(j) + y2(k) + z2(l) < 1.0D0) THEN
          prim2(irho2,j,k,l) = 1.0D0
          prim2(itau2,j,k,l) = 0.1795D0
        ELSEIF(x2(j) + y2(k) + z2(l) > 1.0D0 .AND. x2(j) + y2(k) + z2(l) < 1.5D0) THEN
          prim2(irho2,j,k,l) = 4.0D0    
          prim2(itau2,j,k,l) = 1.0D0
        ELSE
          prim2(irho2,j,k,l) = 1.0D0
          prim2(itau2,j,k,l) = 0.1795D0
        END IF
        prim2(ivel2_x,j,k,l) = 0.0D0
        prim2(ivel2_y,j,k,l) = 0.0D0
        epsilon2(j,k,l) = prim2(itau2,j,k,l) / prim2(irho2,j,k,l) / (ggas2 - 1.0D0)
      END DO
    END DO
  END DO

ENDIF

! set boundary conditions !
call BOUNDARY1D_NM (epsilon2,part,even, even, even, even, even, even)
call BOUNDARYP_NM

! set atmospheric primitive variables !
prim2_a(:) = 0.0D0
eps2_a = 0.0D0

END SUBROUTINE 