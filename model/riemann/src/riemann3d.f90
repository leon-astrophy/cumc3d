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

! 3D Explosion with gamma = 1.4, domain is [0,1] for x, y, and z
IF(test_model == 1) THEN 
  
  ggas = 1.4D0
  total_time = 1.0d-30
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/1,2,1,2,1,2/)

  DO l = 1, nz
    DO k = 1, ny
      DO j = 1, nx
        IF(DSQRT(x(j)**2 + y(k)**2 + z(l)**2) < 0.4D0) THEN
          prim(irho,j,k,l) = 1.0D0       
          prim(itau,j,k,l) = 1.0D0
        ELSE
          prim(irho,j,k,l) = 0.125D0
          prim(itau,j,k,l) = 0.1D0
        END IF
        prim(ivx,j,k,l) = 0.0D0
        prim(ivy,j,k,l) = 0.0D0
        epsilon(j,k,l) = prim(itau,j,k,l) / prim(irho,j,k,l) / (ggas - 1.0D0)
      END DO
    END DO
  END DO

! Oblique Sod shock tube problem
ELSEIF(test_model == 2) THEN

  ggas = (5.0D0/3.0D0)

  DO l = 1, nz
    DO k = 1, ny
      DO j = 1, nx
        IF(x(j) + y(k) + z(l) < 0.5D0) THEN
          prim(irho,j,k,l) = 4.0D0    
          prim(itau,j,k,l) = 1.0D0
        ELSEIF(x(j) + y(k) + z(l) > 0.5D0 .AND. x(j) + y(k) + z(l) < 1.0D0) THEN
          prim(irho,j,k,l) = 1.0D0
          prim(itau,j,k,l) = 0.1795D0
        ELSEIF(x(j) + y(k) + z(l) > 1.0D0 .AND. x(j) + y(k) + z(l) < 1.5D0) THEN
          prim(irho,j,k,l) = 4.0D0    
          prim(itau,j,k,l) = 1.0D0
        ELSE
          prim(irho,j,k,l) = 1.0D0
          prim(itau,j,k,l) = 0.1795D0
        END IF
        prim(ivx,j,k,l) = 0.0D0
        prim(ivy,j,k,l) = 0.0D0
        epsilon(j,k,l) = prim(itau,j,k,l) / prim(irho,j,k,l) / (ggas - 1.0D0)
      END DO
    END DO
  END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

else
	stop "no such test model"
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE 