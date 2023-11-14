!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! 
! This subroutine provides a quick build of initial profiles 
! based on some well known tests 
! It includes: 
! 1. Gresho problem
! 2. Toro blast test
! 3. Toro explosion test
! 4. Implosion test
! 5. Kelvin-Helmholtz Instability
! 6. Rayleigh-Taylor Instability
! 7. 2D Riemann problem case 3
! 8. 2D Riemann problem case 4
! 9. 2D Riemann problem case 6
! 10. 2D Riemann problem case 12
! 11. 2D Riemann problem case 15
! 12. 2D Riemann problem case 17
! 13. MHD Current sheet test
! 14. MHD rotor
! 15. MHD Orszag-Tang Vortex
! 16. MHD blast test
! 17. MHD Kelvin-Helmholtz Instability
! 18. MHD Rayleigh-Taylor Instability
! Written by Leung Shing Chi in 2015 
! Updated by Leung Shing Chi in 2017 
! Included MHD test by Leon Chan in 2022 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Riemann_2d
USE DEFINITION
IMPLICIT NONE
INCLUDE "param.h"

! Integer and real numbers !
INTEGER :: i, j, k, l
REAL*8 :: dummy, r0, r1, r, fr

! Vector potential !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: A_corner

! Magnetic fields (surface) !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: bx_face, by_face

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1. Gresho test with free boundaries everywhere
IF(test_model == 1) THEN

  ggas = 1.4D0
  total_time = 3.0D0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/1,1,1,1,1,1/)

  DO k = 1, ny
    DO j = 1, nx
      dummy = DSQRT((x(j) - 0.5D0)** 2 + (y(k) - 0.5D0)** 2)
      IF(dummy <= 0.2D0) THEN
        prim(itau,j,k,:) = 5.0D0 + 12.5D0 * (dummy ** 2)
        prim(ivx,j,k,:) = -5.0D0 * (y(k) - 0.5D0)
        prim(ivy,j,k,:) = 5.0D0 * (x(j) - 0.5D0)
      ELSEIF(dummy > 0.2D0 .and. dummy <= 0.4D0) THEN
        prim(itau,j,k,:) = 9.0D0 - 4.0D0 * LOG(0.2D0) + 12.5D0 * (dummy ** 2) - 20.0D0 * dummy + 4.0D0 * LOG(dummy)
        prim(ivx,j,k,:) = -(2.0D0 - 5.0D0 * dummy) * (y(k) - 0.5D0) / dummy
        prim(ivy,j,k,:) = (2.0D0 - 5.0D0 * dummy) * (x(j) - 0.5D0) / dummy
      ELSE
        prim(itau,j,k,:) = 3.0D0 + 4.0D0 * LOG(2.0D0)
        prim(ivx,j,k,:) = 0.0D0
        prim(ivy,j,k,:) = 0.0D0
      ENDIF
      prim(irho,j,k,:) = 1.0D0
      epsilon(j,k,:) = prim(itau,j,k,:) / prim(irho,j,k,:) / (ggas - 1.0D0)
    ENDDO
  ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2. Blast test (Toro1997) with gamma = 1.4, x-box = 1, y-box = 1.5
ELSEIF(test_model == 2) THEN

  ggas = 1.4D0
  total_time = 1.0D0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/2,2,2,2,1,1/)

  DO k = 1, ny
    DO j = 1, nx
      dummy = DSQRT((x(j) - 0.5D0)** 2 + (y(k) - 0.5D0)** 2)
      IF(dummy <= 0.1D0) THEN       
        prim(itau,j,k,:) = 10.0D0
      ELSE
        prim(itau,j,k,:) = 1.0D0
      ENDIF
      prim(irho,j,k,:) = 1.0D0
      prim(ivx,j,k,:) = 0.0D0
      prim(ivy,j,k,:) = 0.0D0
      epsilon(j,k,:) = prim(itau,j,k,:) / prim(irho,j,k,:) / (ggas - 1.0D0)
    ENDDO
  ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3. 2D Explosion (Toro1997) with gamma = 1.4, x-box = 1.5, y-box = 1.5
ELSEIF(test_model == 3) THEN

  ggas = 1.4D0
  total_time = 2.667D0
  output_profiletime = total_time/50.0D0

  !IF(coordinate_flag == 0) THEN
    !boundary_flag = (/2,1,2,1,1,1/)
  !ELSEIF(coordinate_flag == 1) THEN
    !boundary_flag = (/3,1,1,1,4,1/)
  !ELSEIF(coordinate_flag == 2) THEN
    !boundary_flag = (/2,1,3,4,1,1/)
  !END IF

  DO k = 1, ny
    DO j = 1, nx
      IF(coordinate_flag == 0) THEN
        dummy = DSQRT(x(j) ** 2 + y(k) ** 2)
      ELSEIF(coordinate_flag == 1) THEN
        dummy = DSQRT(x(j) ** 2 + z(l) ** 2)
      ELSEIF(coordinate_flag == 2) THEN
        dummy = x(j)
      END IF
      IF(dummy < 0.4D0) THEN
        prim(irho,j,k,:) = 1.0D0       
        prim(itau,j,k,:) = 1.0D0
      ELSE
        prim(irho,j,k,:) = 0.125D0
        prim(itau,j,k,:) = 0.1D0
      ENDIF
      prim(ivx,j,k,:) = 0.0D0
      prim(ivy,j,k,:) = 0.0D0
      epsilon(j,k,:) = prim(itau,j,k,:) / prim(irho,j,k,:) / (ggas - 1.0D0)
    ENDDO
  ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 4. Implosion test, http://www-troja.fjfi.cvut.cz/~liska/CompareEuler/compare8/node38_mn.html 
! x-box = 0.3, y-box = 0.3
ELSEIF(test_model == 4) THEN

  ggas = 1.4D0
  total_time = 2.5D0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/2,2,2,2,1,1/)

  DO k = 1, ny
    DO j = 1, nx
      dummy = x(j) + y(k)
			IF (dummy > 0.15) THEN
				prim(irho,j,k,:) = 1.0D0
				prim(itau,j,k,:) = 1.0D0
			ELSE
				prim(irho,j,k,:) = 0.125D0
				prim(itau,j,k,:) = 0.14D0
			END IF
			epsilon(j,k,:) = prim(itau,j,k,:) / prim(irho,j,k,:) / (ggas - 1.0D0)
			prim(ivx,j,k,:) = 0.0D0
			prim(ivy,j,k,:) = 0.0D0
		END DO
	END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Kelvin-Helmholtz test
elseif(test_model == 5) then

  ggas = 1.4D0
  total_time = 5.0d0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/0,0,0,0,1,1/)

  DO k = 1, ny
    DO j = 1, nx
      if(y(k) > 0.75 .OR. y(k) < 0.25) then
        prim(irho,j,k,:) = 1.0D0 
        call random_number(dummy)   
        prim(ivx,j,k,:) = -0.5D0 + 0.01D0*(dummy - 0.5D0)
      else        
        prim(irho,j,k,:) = 2.0D0
        call random_number(dummy) 
        prim(ivx,j,k,:) = 0.5D0 + 0.01D0*(dummy - 0.5D0)
      endif
      call random_number(dummy)
      prim(ivy,j,k,:) = 0.01D0*(dummy - 0.5D0)
      prim(itau,j,k,:) = 2.5D0
      epsilon(j,k,:) = prim(itau,j,k,:) / prim(irho,j,k,:) / (ggas - 1.0D0)
    enddo
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Rayleigh-Taylor, x-box = 0.5, y-box = 1.5
elseif(test_model == 6) then

  ggas = 1.4D0
  total_time = 12.75d0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/0,0,2,2,1,1/)

  DO k = 1, ny
    DO j = 1, nx
      IF(y(k) < 0.75d0) THEN
        prim(irho,j,k,:) = 1.0D0
      ELSE
        prim(irho,j,k,:) = 2.0D0
      ENDIF
      prim(itau,j,k,:) = 2.5D0 - 0.1D0 * prim(irho,j,k,:) * y(k)
      prim(ivx,j,k,:) = 0.0D0
      prim(ivy,j,k,:) = 0.025D0 * (1.0D0 + DCOS(4.0D0 * pi * (x(j) - 0.25D0))) * (1.0D0 + DCOS(3.0D0 * pi * (y(k) - 0.75D0)))
      epsilon(j,k,:) = prim(itau,j,k,:) / prim(irho,j,k,:) / (ggas - 1.0D0)
    ENDDO
  ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2D Riemann Problem Test 3 with gamma = 1.4
ELSEIF(test_model == 7) THEN

  ggas = 1.4D0
  total_time = 0.3D0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/1,1,1,1,1,1/)

  DO k = 1, ny
    DO j = 1, nx
      IF(x(j) > 0.0D0 .AND. x(j) < 0.5D0 .AND. y(k) > 0.0D0 .AND. y(k) < 0.5D0) THEN
        prim(irho,j,k,:) = 0.138D0       
        prim(itau,j,k,:) = 0.0290D0
        prim(ivx,j,k,:) = 1.2060D0
        prim(ivy,j,k,:) = 1.2060D0
      ELSEIF(x(j) > 0.5D0 .AND. x(j) < 1.0D0 .AND. y(k) > 0.0D0 .AND. y(k) < 0.5D0) THEN
        prim(irho,j,k,:) = 0.5323D0       
        prim(itau,j,k,:) = 0.3D0
        prim(ivx,j,k,:) = 0.0D0
        prim(ivy,j,k,:) = 1.2060D0
      ELSEIF(x(j) > 0.0D0 .AND. x(j) < 0.5D0 .AND. y(k) > 0.5D0 .AND. y(k) < 1.0D0) THEN
        prim(irho,j,k,:) = 0.5323D0       
        prim(itau,j,k,:) = 0.3D0
        prim(ivx,j,k,:) = 1.2060D0
        prim(ivy,j,k,:) = 0.0D0
      ELSEIF(x(j) > 0.5D0 .AND. x(j) < 1.0D0 .AND. y(k) > 0.5D0 .AND. y(k) < 1.0D0) THEN
        prim(irho,j,k,:) = 1.5D0  
        prim(itau,j,k,:) = 1.5D0 
        prim(ivx,j,k,:) = 0.0D0
        prim(ivy,j,k,:) = 0.0D0
      END IF
      epsilon(j,k,:) = prim(itau,j,k,:) / prim(irho,j,k,:) / (ggas - 1.0D0)
    ENDDO
  ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2D Riemann Problem Test 4 with gamma = 1.4
ELSEIF(test_model == 8) THEN

  ggas = 1.4D0
  total_time = 0.25D0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/1,1,1,1,1,1/)

  DO k = 1, ny
    DO j = 1, nx
      IF(x(j) > 0.0D0 .AND. x(j) < 0.5D0 .AND. y(k) > 0.0D0 .AND. y(k) < 0.5D0) THEN
        prim(itau,j,k,:) = 1.1D0
        prim(irho,j,k,:) = 1.1D0       
        prim(ivx,j,k,:) = 0.8939D0
        prim(ivy,j,k,:) = 0.8939D0
      ELSEIF(x(j) > 0.5D0 .AND. x(j) < 1.0D0 .AND. y(k) > 0.0D0 .AND. y(k) < 0.5D0) THEN
        prim(itau,j,k,:) = 0.35D0
        prim(irho,j,k,:) = 0.5065D0       
        prim(ivx,j,k,:) = 0.0D0
        prim(ivy,j,k,:) = 0.8939D0
      ELSEIF(x(j) > 0.0D0 .AND. x(j) < 0.5D0 .AND. y(k) > 0.5D0 .AND. y(k) < 1.0D0) THEN
        prim(itau,j,k,:) = 0.35D0
        prim(irho,j,k,:) = 0.5065D0       
        prim(ivx,j,k,:) = 0.8939D0
        prim(ivy,j,k,:) = 0.0D0
      ELSEIF(x(j) > 0.5D0 .AND. x(j) < 1.0D0 .AND. y(k) > 0.5D0 .AND. y(k) < 1.0D0) THEN
        prim(itau,j,k,:) = 1.1D0  
        prim(irho,j,k,:) = 1.1D0 
        prim(ivx,j,k,:) = 0.0D0
        prim(ivy,j,k,:) = 0.0D0
      END IF
      epsilon(j,k,:) = prim(itau,j,k,:) / prim(irho,j,k,:) / (ggas - 1.0D0)
    ENDDO
  ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2D Riemann Problem Test 6 with gamma = 1.4
ELSEIF(test_model == 9) THEN

  ggas = 1.4D0
  total_time = 0.3D0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/1,1,1,1,1,1/)

  DO k = 1, ny
    DO j = 1, nx
      IF(x(j) > 0.0D0 .AND. x(j) < 0.5D0 .AND. y(k) > 0.0D0 .AND. y(k) < 0.5D0) THEN
        prim(itau,j,k,:) = 1.0D0       
        prim(irho,j,k,:) = 1.0D0
        prim(ivx,j,k,:) = -0.75D0
        prim(ivy,j,k,:) = 0.5D0
      ELSEIF(x(j) > 0.5D0 .AND. x(j) < 1.0D0 .AND. y(k) > 0.0D0 .AND. y(k) < 0.5D0) THEN
        prim(itau,j,k,:) = 1.0D0       
        prim(irho,j,k,:) = 3.0D0
        prim(ivx,j,k,:) = -0.75D0
        prim(ivy,j,k,:) = -0.5D0
      ELSEIF(x(j) > 0.0D0 .AND. x(j) < 0.5D0 .AND. y(k) > 0.5D0 .AND. y(k) < 1.0D0) THEN
        prim(itau,j,k,:) = 1.0D0       
        prim(irho,j,k,:) = 2.0D0
        prim(ivx,j,k,:) = 0.75D0
        prim(ivy,j,k,:) = 0.5D0
      ELSEIF(x(j) > 0.5D0 .AND. x(j) < 1.0D0 .AND. y(k) > 0.5D0 .AND. y(k) < 1.0D0) THEN
        prim(itau,j,k,:) = 1.0D0  
        prim(irho,j,k,:) = 1.0D0 
        prim(ivx,j,k,:) = 0.75D0
        prim(ivy,j,k,:) = -0.5D0
      END IF
      epsilon(j,k,:) = prim(itau,j,k,:) / prim(irho,j,k,:) / (ggas - 1.0D0)
    ENDDO
  ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2D Riemann Problem Test 12 with gamma = 1.4
ELSEIF(test_model == 10) THEN

  ggas = 1.4D0
  total_time = 0.25D0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/1,1,1,1,1,1/)

  DO k = 1, ny
    DO j = 1, nx
      IF(x(j) > 0.0D0 .AND. x(j) < 0.5D0 .AND. y(k) > 0.0D0 .AND. y(k) < 0.5D0) THEN
        prim(itau,j,k,:) = 1.0D0       
        prim(irho,j,k,:) = 0.8D0
        prim(ivx,j,k,:) = 0.0D0
        prim(ivy,j,k,:) = 0.0D0
      ELSEIF(x(j) > 0.5D0 .AND. x(j) < 1.0D0 .AND. y(k) > 0.0D0 .AND. y(k) < 0.5D0) THEN
        prim(itau,j,k,:) = 1.0D0       
        prim(irho,j,k,:) = 1.0D0
        prim(ivx,j,k,:) = 0.0D0
        prim(ivy,j,k,:) = 0.7276D0
      ELSEIF(x(j) > 0.0D0 .AND. x(j) < 0.5D0 .AND. y(k) > 0.5D0 .AND. y(k) < 1.0D0) THEN
        prim(itau,j,k,:) = 1.0D0       
        prim(irho,j,k,:) = 1.0D0
        prim(ivx,j,k,:) = 0.7276D0
        prim(ivy,j,k,:) = 0.5D0
      ELSEIF(x(j) > 0.5D0 .AND. x(j) < 1.0D0 .AND. y(k) > 0.5D0 .AND. y(k) < 1.0D0) THEN
        prim(itau,j,k,:) = 0.4D0  
        prim(irho,j,k,:) = 0.5313D0 
        prim(ivx,j,k,:) = 0.0D0
        prim(ivy,j,k,:) = 0.0D0
      END IF
      epsilon(j,k,:) = prim(itau,j,k,:) / prim(irho,j,k,:) / (ggas - 1.0D0)
    ENDDO
  ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2D Riemann Problem Test 15 with gamma = 1.4
ELSEIF(test_model == 11) THEN

  ggas = 1.4D0
  total_time = 0.2D0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/1,1,1,1,1,1/)

  DO k = 1, ny
    DO j = 1, nx
      IF(x(j) > 0.0D0 .AND. x(j) < 0.5D0 .AND. y(k) > 0.0D0 .AND. y(k) < 0.5D0) THEN
        prim(itau,j,k,:) = 0.4D0       
        prim(irho,j,k,:) = 0.8D0
        prim(ivx,j,k,:) = 0.1D0
        prim(ivy,j,k,:) = -0.3D0
      ELSEIF(x(j) > 0.5D0 .AND. x(j) < 1.0D0 .AND. y(k) > 0.0D0 .AND. y(k) < 0.5D0) THEN
        prim(itau,j,k,:) = 0.4D0       
        prim(irho,j,k,:) = 0.5313D0
        prim(ivx,j,k,:) = 0.1D0
        prim(ivy,j,k,:) = 0.4276D0
      ELSEIF(x(j) > 0.0D0 .AND. x(j) < 0.5D0 .AND. y(k) > 0.5D0 .AND. y(k) < 1.0D0) THEN
        prim(itau,j,k,:) = 0.4D0       
        prim(irho,j,k,:) = 0.5197D0
        prim(ivx,j,k,:) = -0.6259D0
        prim(ivy,j,k,:) = -0.3D0
      ELSEIF(x(j) > 0.5D0 .AND. x(j) < 1.0D0 .AND. y(k) > 0.5D0 .AND. y(k) < 1.0D0) THEN
        prim(itau,j,k,:) = 1.0D0  
        prim(irho,j,k,:) = 1.0D0 
        prim(ivx,j,k,:) = 0.1D0
        prim(ivy,j,k,:) = -0.3D0
      END IF
      epsilon(j,k,:) = prim(itau,j,k,:) / prim(irho,j,k,:) / (ggas - 1.0D0)
    ENDDO
  ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2D Riemann Problem Test 17 with gamma = 1.4
ELSEIF(test_model == 12) THEN

  ggas = 1.4D0
  total_time = 0.3D0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/1,1,1,1,1,1/)

  DO k = 1, ny
    DO j = 1, nx
      IF(x(j) > 0.0D0 .AND. x(j) < 0.5D0 .AND. y(k) > 0.0D0 .AND. y(k) < 0.5D0) THEN
        prim(itau,j,k,:) = 0.4D0       
        prim(irho,j,k,:) = 1.0625D0
        prim(ivx,j,k,:) = 0.0D0
        prim(ivy,j,k,:) = 0.2145D0
      ELSEIF(x(j) > 0.5D0 .AND. x(j) < 1.0D0 .AND. y(k) > 0.0D0 .AND. y(k) < 0.5D0) THEN
        prim(itau,j,k,:) = 0.4D0       
        prim(irho,j,k,:) = 0.5197D0
        prim(ivx,j,k,:) = 0.0D0
        prim(ivy,j,k,:) = -1.1259D0
      ELSEIF(x(j) > 0.0D0 .AND. x(j) < 0.5D0 .AND. y(k) > 0.5D0 .AND. y(k) < 1.0D0) THEN
        prim(itau,j,k,:) = 1.0D0       
        prim(irho,j,k,:) = 2.0D0
        prim(ivx,j,k,:) = 0.0D0
        prim(ivy,j,k,:) = -0.3D0
      ELSEIF(x(j) > 0.5D0 .AND. x(j) < 1.0D0 .AND. y(k) > 0.5D0 .AND. y(k) < 1.0D0) THEN
        prim(itau,j,k,:) = 1.0D0  
        prim(irho,j,k,:) = 1.0D0 
        prim(ivx,j,k,:) = 0.0D0
        prim(ivy,j,k,:) = -0.4D0
      END IF
      epsilon(j,k,:) = prim(itau,j,k,:) / prim(irho,j,k,:) / (ggas - 1.0D0)
    ENDDO
  ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MHD current sheet test
ELSEIF(test_model == 13) THEN

  ggas = 5.0d0/3.0d0
  total_time = 10.0D0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/0,0,0,0,1,1/)

  DO k = 0, ny
    DO j = 0, nx
      IF(x(j) <= 0.25d0 .OR. x(j) >= 0.75d0) THEN
        prim(iby,j,k,:) = 1.0D0/DSQRT(4.0D0*pi)
      ELSE
        prim(iby,j,k,:) = -1.0D0/DSQRT(4.0D0*pi)
      END IF
      prim(irho,j,k,:) = 1.0D0
      prim(itau,j,k,:) = 0.1D0*0.5D0
      prim(ivx,j,k,:) = 0.1D0*DSIN(2.0D0*pi*(y(k) - 0.5D0))
      epsilon(j,k,:) = prim(itau,j,k,:) / prim(irho,j,k,:) / (ggas - 1.0D0)
    ENDDO
  ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MHD rotor test
ELSEIF(test_model == 14) THEN

  ggas = (5.0D0/3.0D0)
  r0 = 0.1D0
  r1 = 0.115D0
  total_time = 0.295D0
  output_profiletime = total_time/10.0D0
  
  !boundary_flag = (/1,1,1,1,1,1/)

  DO k = 0, ny
    DO j = 0, nx
      r = DSQRT((x(j) - 0.5D0)**2 + (y(k) - 0.5D0)**2)
      fr = (r1 - r)/(r1 - r0)
      If(r <= r0) THEN
        prim(irho,j,k,:) = 10.0D0
        prim(ivx,j,k,:) = -1.0D0*(y(k) - 0.5D0)/r0
        prim(ivy,j,k,:) = 1.0D0*(x(j) - 0.5D0)/r0
      ELSEIF(r >= r1) THEN
        prim(irho,j,k,:) = 1.0D0
        prim(ivx,j,k,:) = 0.0D0
        prim(ivy,j,k,:) = 0.0D0
      ELSE
        prim(irho,j,k,:) = 1.0D0 + 9.0d0*fr
        prim(ivx,j,k,:) = -fr*1.0D0*(y(k) - 0.5D0)/r
        prim(ivy,j,k,:) = fr*1.0D0*(x(j) - 0.5D0)/r
      END IF
      prim(itau,j,k,:) = 0.5D0
      prim(ibx,j,k,:) = 2.5D0/DSQRT(4.0D0*pi)
      epsilon(j,k,:) = prim(itau,j,k,:) / prim(irho,j,k,:) / (ggas - 1.0D0)
    ENDDO
  ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MHD Orszag-Tang Vortex test
ELSEIF(test_model == 15) THEN

  ggas = (5.0D0/3.0D0)
  total_time = 1.0D0
  output_profiletime = total_time/50.0D0

  !boundary_flag = (/0,0,0,0,1,1/)

  DO k = 0, ny
    DO j = 0, nx
      prim(irho,j,k,:) = 25.0D0/36.0D0/pi
      prim(itau,j,k,:) = 5.0D0/12.0D0/pi
      prim(ivx,j,k,:) = - DSIN(2.0D0*pi*y(k))
      prim(ivy,j,k,:) = DSIN(2.0D0*pi*x(j))
      prim(ibx,j,k,:) = - DSIN(2.0D0*pi*y(k))/DSQRT(4.0D0*pi)
      prim(iby,j,k,:) = DSIN(4.0D0*pi*x(j))/DSQRT(4.0D0*pi)
      epsilon(j,k,:) = prim(itau,j,k,:) / prim(irho,j,k,:) / (ggas - 1.0D0)
    ENDDO
  ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MHD blast test
ELSEIF(test_model == 16) THEN

  ggas = (5.0D0/3.0D0)
  total_time = 0.2D0
  output_profiletime = total_time/50.0D0

  !boundary_flag = (/0,0,0,0,1,1/)

  DO k = 0, ny
    DO j = 0, nx
      dummy = DSQRT((x(j) - 0.5D0)**2 + (y(k) - 0.5D0)**2)
      IF(dummy <= 0.1D0) THEN
        prim(itau,j,k,:) = 10.0D0
      ELSE
        prim(itau,j,k,:) = 0.1D0
      END If
      prim(irho,j,k,:) = 1.0D0
      prim(ibx,j,k,:) = 1.0D0/DSQRT(2)
      prim(iby,j,k,:) = 1.0D0/DSQRT(2)
      epsilon(j,k,:) = prim(itau,j,k,:) / prim(irho,j,k,:) / (ggas - 1.0D0)
    ENDDO
  ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MHD Kelvin-Helmholtz test
elseif(test_model == 17) then

  ggas = 1.4D0
  total_time = 5.0d0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/0,0,0,0,1,1/)

  DO k = 0, ny
    DO j = 0, nx
      if(y(k) > 0.75 .OR. y(k) < 0.25) then
        prim(irho,j,k,:) = 1.0D0    
        call random_number(dummy)   
        prim(ivx,j,k,:) = 0.5D0 + 0.01D0*dummy
      else        
        prim(irho,j,k,:) = 2.0D0
        call random_number(dummy) 
        prim(ivx,j,k,:) = -0.5D0 + 0.01D0*dummy
      endif
      prim(itau,j,k,:) = 2.5D0
      call random_number(dummy)
      prim(ivy,j,k,:) = 0.01D0*dummy
      epsilon(j,k,:) = prim(itau,j,k,:) / prim(irho,j,k,:) / (ggas - 1.0D0)
      prim(ibx,j,k,:) = 0.2d0
    enddo
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MHD Rayleigh-Taylor, x-box = 1.0, y-box = 1.0
elseif(test_model == 18) then

  ggas = 1.4D0
  total_time = 12.5D0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/0,0,2,2,1,1/)

  DO k = 0, ny
    DO j = 0, nx
      IF(y(k) <= 0.5d0) THEN
        prim(irho,j,k,:) = 1.0D0
      ELSE
        prim(irho,j,k,:) = 2.0D0
      ENDIF
      prim(itau,j,k,:) = 2.5D0 - 0.1D0 * prim(irho,j,k,:) * (y(k) - 0.5D0)
      prim(ivx,j,k,:) = 0.0D0
      call random_number(dummy)
      prim(ivy,j,k,:) = 0.01D0*dummy*(1.0d0 + DCOS(8.0d0*pi*(y(k) - 0.5D0)/3.0d0))/2.0d0
      epsilon(j,k,:) = prim(itau,j,k,:) / prim(irho,j,k,:) / (ggas - 1.0D0)
      prim(ibx,j,k,:) = 0.0125D0
    ENDDO
  ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MHD advection of current loop , -1 < x < 1, -0.5 < y < 0.5 !
elseif(test_model == 19) then

  ggas = 5.0d0/3.0d0
  total_time = 2.0d0
  output_profiletime = total_time/50.0D0

  !boundary_flag = (/0,0,0,0,1,1/)

  ! Allocate ! 
  ALLOCATE(A_corner(-2:nx+3,-2:ny+3,-2:nz+3))

  A_corner = 0.0d0
  DO k = 0, ny
    DO j = 0, nx
      prim(irho,j,k,:) = 1.0d0
      prim(itau,j,k,:) = 1.0d0
      prim(ivx,j,k,:) = 2.0d0
      prim(ivy,j,k,:) = 1.0d0
      epsilon(j,k,:) = prim(itau,j,k,:) / prim(irho,j,k,:) / (ggas - 1.0D0)

      ! cell corner coordinate !
      r = DSQRT(xF(j)**2 + yF(k)**2)
      A_corner(j,k,:) = 0.001D0*MAX(0.0d0, (0.3d0 - r))
    ENDDO
  ENDDO

  ! Face-centered magnetic field !
  DO k = 0, ny
    DO j = 0, nx
      prim(ibx,j,k,:) = (A_corner(j,k,:) - A_corner(j,k-1,:))/dy(k)
      prim(iby,j,k,:) = -(A_corner(j,k,:) - A_corner(j-1,k,:))/dx(j)
    ENDDO
  ENDDO

  ! Deallocate ! 
  DEALLOCATE(A_corner)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Error !
  
else
	stop "no such test model"
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE 
