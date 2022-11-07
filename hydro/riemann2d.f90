!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! 
! This subroutine provides a quick build of initial profiles 
! based on some well known tests 
! It includes: 
! 1. Sedov Spherical test
! 2. 2D explosion
! 3. Kelvin-helmholtz
! 4. Rayleigh-Taylor
! 5. Kelvin-helmholtz, another setting
! 6. Gresho problem
! 7. Implosion
! 8 - 13. 2D Riemann problem
! Written by Leung Shing Chi in 2015 
! Updated by Leung Shing Chi in 2017 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Riemann_2d
USE DEFINITION
IMPLICIT NONE

! Integer and real numbers !
INTEGER :: i, j, k
REAL (DP) :: dummy

! 2D shock test with gamma = 1.4
IF(test_model == 1) THEN

   ggas2 = 1.4E0_DP

   DO k = 1, ny_2
      DO j = 1, nx_2
         if(DSQRT((x2(j) - 0.5D0) ** 2 + (y2(k) - 0.75D0) ** 2) < 0.1D0) then       
            prim2(irho2,j,k,:) = 1.0D0       
            prim2(itau2,j,k,:) = 10.0D0
            prim2(ivel2_x,j,k,:) = 0.0D0
            prim2(ivel2_y,j,k,:) = 0.0D0
            epsilon2(j,k,:) = prim2(itau2,j,k,:) / prim2(irho2,j,k,:) / (ggas2 - 1.0D0)
         else
            prim2(irho2,j,k,:) = 1.0D0
            prim2(itau2,j,k,:) = 1.0D0
            prim2(ivel2_x,j,k,:) = 0.0D0
            prim2(ivel2_y,j,k,:) = 0.0D0
            epsilon2(j,k,:) = prim2(itau2,j,k,:) / prim2(irho2,j,k,:) / (ggas2 - 1.0D0)
         endif
      ENDDO
   ENDDO

! 2D Explosion (Toro1997) with gamma = 1.4
ELSEIF(test_model == 2) THEN

   ggas2 = 1.4E0_DP

   DO k = 1, ny_2
      DO j = 1, nx_2
         IF(DSQRT(x2(j) ** 2 + y2(k) ** 2) < 0.4D0) THEN
            prim2(irho2,j,k,:) = 1.0D0       
            prim2(itau2,j,k,:) = 1.0D0
            prim2(ivel2_x,j,k,:) = 0.0D0
            prim2(ivel2_y,j,k,:) = 0.0D0
            epsilon2(j,k,:) = prim2(itau2,j,k,:) / prim2(irho2,j,k,:) / (ggas2 - 1.0D0)
         ELSE
            prim2(irho2,j,k,:) = 0.125D0
            prim2(itau2,j,k,:) = 0.1D0
            prim2(ivel2_x,j,k,:) = 0.0D0
            prim2(ivel2_y,j,k,:) = 0.0D0
            epsilon2(j,k,:) = prim2(itau2,j,k,:) / prim2(irho2,j,k,:) / (ggas2 - 1.0D0)
         ENDIF
      ENDDO
   ENDDO

! Kelvin-Helmholtz with gamma = 1.4
elseif(test_model == 3) then

   ggas2 = 1.4E0_DP

   DO k = 1, ny_2
      DO j = 1, nx_2
         if(k <= ny_2 * 3 / 4 .and. k >= ny_2 / 4) then
            prim2(irho2,j,k,:) = 2.0D0 
            prim2(itau2,j,k,:) = 2.5D0      
            prim2(ivel2_x,j,k,:) = -0.5D0 
            prim2(ivel2_y,j,k,:) = 0.05D0 * DSIN(4.0D0 * pi * x2(j))
            epsilon2(j,k,:) = prim2(itau2,j,k,:) / prim2(irho2,j,k,:) / (ggas2 - 1.0D0)
         else                  
            prim2(irho2,j,k,:) = 1.0D0
            prim2(itau2,j,k,:) = 2.5D0
            prim2(ivel2_x,j,k,:) = 0.5D0
            prim2(ivel2_y,j,k,:) = 0.05D0 * DSIN(4.0D0 * pi * x2(j))
            epsilon2(j,k,:) = prim2(itau2,j,k,:) / prim2(irho2,j,k,:) / (ggas2 - 1.0D0)
         endif
      enddo
   enddo

! Rayleigh-Taylor
ELSEIF(test_model == 4) THEN

   ggas2 = 1.4E0_DP

   DO k = 1, ny_2
      DO j = 1, nx_2
         IF(k <= ny_2 * 1 / 2) THEN
            prim2(irho2,j,k,:) = 1.0D0
            prim2(itau2,j,k,:) = 2.5D0 - 0.1D0 * prim2(irho2,j,k,:) * y2(k)
            prim2(ivel2_x,j,k,:) = 0.0D0
            prim2(ivel2_y,j,k,:) = 0.025D0 * (1.0D0 + DCOS(4.0D0 * pi * (x2(j) - 0.25D0))) * &
                                    (1.0D0 + DCOS(3.0D0 * pi * (y2(k) - 0.75D0)))
            epsilon2(j,k,:) = prim2(itau2,j,k,:) / prim2(irho2,j,k,:) / (ggas2 - 1.0D0)
         ELSE
            prim2(irho2,j,k,:) = 2.0D0
            prim2(itau2,j,k,:) = 2.5D0 - 0.1D0 * prim2(irho2,j,k,:) * y2(k)
            prim2(ivel2_x,j,k,:) = 0.0D0
            prim2(ivel2_y,j,k,:) = 0.025D0 * (1.0D0 + DCOS(4.0D0 * pi * (x2(j) - 0.25D0))) * &
                                    (1.0D0 + DCOS(3.0D0 * pi * (y2(k) - 0.75D0)))
            epsilon2(j,k,:) = prim2(itau2,j,k,:) / prim2(irho2,j,k,:) / (ggas2 - 1.0D0)
         ENDIF
      ENDDO
   ENDDO

! Kelvin-Helmholtz, another setting 
ELSEIF(test_model == 5) THEN

   ! A unit box is assumed
   ggas2 = 1.4E0_DP

   DO k = 1, ny_2
      DO j = 1, nx_2
         IF(k <= ny_2 / 4) THEN
   	      prim2(irho2,j,k,:) = 1.0D0 + 0.5D0 * EXP((y2(k) - 0.25D0) / 0.025D0)
	         prim2(itau2,j,k,:) = 2.5D0
	         prim2(ivel2_x,j,k,:) = 0.5D0 - 0.5D0 * EXP((y2(k) - 0.25D0) / 0.025D0)
	         prim2(ivel2_y,j,k,:) = 0.01D0 * SIN(4.0D0 * pi * x2(j))
	         epsilon2(j,k,:) = prim2(itau2,j,k,:) / prim2(irho2,j,k,:) / (ggas2 - 1.0D0)
	      ELSEIF(k > ny_2 / 4 .and. k <= ny_2 / 2) THEN
	         prim2(irho2,j,k,:) = 2.0D0 - 0.5D0 * EXP((-y2(k) + 0.25D0) / 0.025D0)
            prim2(itau2,j,k,:) = 2.5D0
            prim2(ivel2_x,j,k,:) = -0.5D0 + 0.5D0 * EXP((-y2(k) + 0.25D0) / 0.025D0)
            prim2(ivel2_y,j,k,:) = 0.01D0 * SIN(4.0D0 * pi * x2(j))
            epsilon2(j,k,:) = prim2(itau2,j,k,:) / prim2(irho2,j,k,:) / (ggas2 - 1.0D0)
	      ELSEIF(k > ny_2 / 2 .and. k <= 3 * ny_2 / 4) THEN
	         prim2(irho2,j,k,:) = 2.0D0 - 0.5D0 * EXP(-(0.75D0 - y2(k)) / 0.025D0)
            prim2(itau2,j,k,:) = 2.5D0
            prim2(ivel2_x,j,k,:) = -0.5D0 + 0.5D0 * EXP((-(0.75D0 - y2(k))) / 0.025D0)
            prim2(ivel2_y,j,k,:) = 0.01D0 * SIN(4.0D0 * pi * x2(j))
            epsilon2(j,k,:) = prim2(itau2,j,k,:) / prim2(irho2,j,k,:) / (ggas2 - 1.0D0)
	      ELSE
	         prim2(irho2,j,k,:) = 1.0D0 + 0.5D0 * EXP(-(y2(k) - 0.75D0) / 0.025D0)
            prim2(itau2,j,k,:) = 2.5D0
            prim2(ivel2_x,j,k,:) = 0.5D0 - 0.5D0 * EXP((-(y2(k) - 0.75D0)) / 0.025D0)
            prim2(ivel2_y,j,k,:) = 0.01D0 * SIN(4.0D0 * pi * x2(j))
            epsilon2(j,k,:) = prim2(itau2,j,k,:) / prim2(irho2,j,k,:) / (ggas2 - 1.0D0)
	      ENDIF
      ENDDO
   ENDDO

! Gresho
ELSEIF(test_model == 6) THEN

   ggas2 = 1.4E0_DP

   DO k = 1, ny_2
      DO j = 1, nx_2

         dummy = DSQRT((x2(j) - 0.5D0)** 2 + (y2(k) - 0.5D0)** 2)

         IF(dummy <= 0.2D0) THEN
            prim2(irho2,j,k,:) = 1.0D0
            prim2(itau2,j,k,:) = 5.0D0 + 12.5D0 * (dummy ** 2)
            prim2(ivel2_x,j,k,:) = -5.0D0 * (y2(k) - 0.5D0)
            prim2(ivel2_y,j,k,:) = 5.0D0 * (x2(j) - 0.5D0)
            epsilon2(j,k,:) = prim2(itau2,j,k,:) / prim2(irho2,j,k,:) / (ggas2 - 1.0D0)
         ELSEIF(dummy > 0.2D0 .and. dummy <= 0.4D0) THEN
            prim2(irho2,j,k,:) = 1.0D0
            prim2(itau2,j,k,:) = 9.0D0 - 4.0D0 * LOG(0.2D0) + 12.5D0 * (dummy ** 2) - 20.0D0 * dummy + 4.0D0 * LOG(dummy)
            prim2(ivel2_x,j,k,:) = -(2.0D0 - 5.0D0 * dummy) * (y2(k) - 0.5D0) / dummy
            prim2(ivel2_y,j,k,:) = (2.0D0 - 5.0D0 * dummy) * (x2(j) - 0.5D0) / dummy
            epsilon2(j,k,:) = prim2(itau2,j,k,:) / prim2(irho2,j,k,:) / (ggas2 - 1.0D0)
         ELSE
            prim2(irho2,j,k,:) = 1.0D0
            prim2(itau2,j,k,:) = 3.0D0 + 4.0D0 * LOG(2.0D0)
            prim2(ivel2_x,j,k,:) = 0.0D0
            prim2(ivel2_y,j,k,:) = 0.0D0
            epsilon2(j,k,:) = prim2(itau2,j,k,:) / prim2(irho2,j,k,:) / (ggas2 - 1.0D0)
         ENDIF

      ENDDO
   ENDDO

! Implosion
ELSEIF(test_model == 7) THEN

	ggas2 = 1.4E0_DP

   DO k = 1, ny_2
      DO j = 1, nx_2
			IF ( x2(j) + y2(k) > 0.15 ) THEN
				prim2(irho2,j,k,:) = 1.0E0_DP
				prim2(itau2,j,k,:) = 1.0E0_DP
			ELSE
				prim2(irho2,j,k,:) = 0.125E0_DP
				prim2(itau2,j,k,:) = 0.14E0_DP
			END IF
			epsilon2(j,k,:) = prim2(itau2,j,k,:) / prim2(irho2,j,k,:) / (ggas2 - 1.0D0)
			prim2(ivel2_x,j,k,:) = 0.0D0
			prim2(ivel2_y,j,k,:) = 0.0D0
		END DO
	END DO

! 2D Riemann Problem Test with gamma = 1.4
ELSEIF(test_model == 8) THEN

   ggas2 = 1.4D0

   DO k = 1, ny_2
      DO j = 1, nx_2
	      IF(x2(j) > 0.0D0 .AND. x2(j) < 0.5D0 .AND. y2(k) > 0.0D0 .AND. y2(k) < 0.5D0) THEN
            prim2(irho2,j,k,:) = 0.138D0       
            prim2(itau2,j,k,:) = 0.0290D0
            prim2(ivel2_x,j,k,:) = 1.2060D0
            prim2(ivel2_y,j,k,:) = 1.2060D0
	      ELSEIF(x2(j) > 0.5D0 .AND. x2(j) < 1.0D0 .AND. y2(k) > 0.0D0 .AND. y2(k) < 0.5D0) THEN
            prim2(irho2,j,k,:) = 0.5323D0       
            prim2(itau2,j,k,:) = 0.3D0
            prim2(ivel2_x,j,k,:) = 0.0D0
            prim2(ivel2_y,j,k,:) = 1.2060D0
	      ELSEIF(x2(j) > 0.0D0 .AND. x2(j) < 0.5D0 .AND. y2(k) > 0.5D0 .AND. y2(k) < 1.0D0) THEN
            prim2(irho2,j,k,:) = 0.5323D0       
            prim2(itau2,j,k,:) = 0.3D0
            prim2(ivel2_x,j,k,:) = 1.2060D0
            prim2(ivel2_y,j,k,:) = 0.0D0
	      ELSEIF(x2(j) > 0.5D0 .AND. x2(j) < 1.0D0 .AND. y2(k) > 0.5D0 .AND. y2(k) < 1.0D0) THEN
            prim2(irho2,j,k,:) = 1.5D0  
            prim2(itau2,j,k,:) = 1.5D0 
            prim2(ivel2_x,j,k,:) = 0.0D0
            prim2(ivel2_y,j,k,:) = 0.0D0
	      END IF
	      epsilon2(j,k,:) = prim2(itau2,j,k,:) / prim2(irho2,j,k,:) / (ggas2 - 1.0D0)
      ENDDO
   ENDDO

! 2D Riemann Problem Test with gamma = 1.4
ELSEIF(test_model == 9) THEN

   ggas2 = 1.4D0

   DO k = 1, ny_2
      DO j = 1, nx_2
	      IF(x2(j) > 0.0D0 .AND. x2(j) < 0.5D0 .AND. y2(k) > 0.0D0 .AND. y2(k) < 0.5D0) THEN
	         prim2(itau2,j,k,:) = 1.1D0
            prim2(irho2,j,k,:) = 1.1D0       
            prim2(ivel2_x,j,k,:) = 0.8939D0
            prim2(ivel2_y,j,k,:) = 0.8939D0
	      ELSEIF(x2(j) > 0.5D0 .AND. x2(j) < 1.0D0 .AND. y2(k) > 0.0D0 .AND. y2(k) < 0.5D0) THEN
 	         prim2(itau2,j,k,:) = 0.35D0
            prim2(irho2,j,k,:) = 0.5065D0       
            prim2(ivel2_x,j,k,:) = 0.0D0
            prim2(ivel2_y,j,k,:) = 0.8939D0
	      ELSEIF(x2(j) > 0.0D0 .AND. x2(j) < 0.5D0 .AND. y2(k) > 0.5D0 .AND. y2(k) < 1.0D0) THEN
	         prim2(itau2,j,k,:) = 0.35D0
            prim2(irho2,j,k,:) = 0.5065D0       
            prim2(ivel2_x,j,k,:) = 0.8939D0
            prim2(ivel2_y,j,k,:) = 0.0D0
	      ELSEIF(x2(j) > 0.5D0 .AND. x2(j) < 1.0D0 .AND. y2(k) > 0.5D0 .AND. y2(k) < 1.0D0) THEN
            prim2(itau2,j,k,:) = 1.1D0  
            prim2(irho2,j,k,:) = 1.1D0 
            prim2(ivel2_x,j,k,:) = 0.0D0
            prim2(ivel2_y,j,k,:) = 0.0D0
	      END IF
	      epsilon2(j,k,:) = prim2(itau2,j,k,:) / prim2(irho2,j,k,:) / (ggas2 - 1.0D0)
      ENDDO
   ENDDO

! 2D Riemann Problem Test with gamma = 1.4
ELSEIF(test_model == 10) THEN

   ggas2 = 1.4D0

   DO k = 1, ny_2
      DO j = 1, nx_2
	      IF(x2(j) > 0.0D0 .AND. x2(j) < 0.5D0 .AND. y2(k) > 0.0D0 .AND. y2(k) < 0.5D0) THEN
            prim2(itau2,j,k,:) = 1.0D0       
            prim2(irho2,j,k,:) = 1.0D0
            prim2(ivel2_x,j,k,:) = -0.75D0
            prim2(ivel2_y,j,k,:) = 0.5D0
	      ELSEIF(x2(j) > 0.5D0 .AND. x2(j) < 1.0D0 .AND. y2(k) > 0.0D0 .AND. y2(k) < 0.5D0) THEN
            prim2(itau2,j,k,:) = 1.0D0       
            prim2(irho2,j,k,:) = 3.0D0
            prim2(ivel2_x,j,k,:) = -0.75D0
            prim2(ivel2_y,j,k,:) = -0.5D0
	      ELSEIF(x2(j) > 0.0D0 .AND. x2(j) < 0.5D0 .AND. y2(k) > 0.5D0 .AND. y2(k) < 1.0D0) THEN
            prim2(itau2,j,k,:) = 1.0D0       
            prim2(irho2,j,k,:) = 2.0D0
            prim2(ivel2_x,j,k,:) = 0.75D0
            prim2(ivel2_y,j,k,:) = 0.5D0
	      ELSEIF(x2(j) > 0.5D0 .AND. x2(j) < 1.0D0 .AND. y2(k) > 0.5D0 .AND. y2(k) < 1.0D0) THEN
            prim2(itau2,j,k,:) = 1.0D0  
            prim2(irho2,j,k,:) = 1.0D0 
            prim2(ivel2_x,j,k,:) = 0.75D0
            prim2(ivel2_y,j,k,:) = -0.5D0
	      END IF
	      epsilon2(j,k,:) = prim2(itau2,j,k,:) / prim2(irho2,j,k,:) / (ggas2 - 1.0D0)
      ENDDO
   ENDDO

! 2D Riemann Problem Test with gamma = 1.4
ELSEIF(test_model == 11) THEN

   ggas2 = 1.4D0

   DO k = 1, ny_2
      DO j = 1, nx_2
	      IF(x2(j) > 0.0D0 .AND. x2(j) < 0.5D0 .AND. y2(k) > 0.0D0 .AND. y2(k) < 0.5D0) THEN
            prim2(itau2,j,k,:) = 1.0D0       
            prim2(irho2,j,k,:) = 0.8D0
            prim2(ivel2_x,j,k,:) = 0.0D0
            prim2(ivel2_y,j,k,:) = 0.0D0
	      ELSEIF(x2(j) > 0.5D0 .AND. x2(j) < 1.0D0 .AND. y2(k) > 0.0D0 .AND. y2(k) < 0.5D0) THEN
            prim2(itau2,j,k,:) = 1.0D0       
            prim2(irho2,j,k,:) = 1.0D0
            prim2(ivel2_x,j,k,:) = 0.0D0
            prim2(ivel2_y,j,k,:) = 0.7276D0
	      ELSEIF(x2(j) > 0.0D0 .AND. x2(j) < 0.5D0 .AND. y2(k) > 0.5D0 .AND. y2(k) < 1.0D0) THEN
            prim2(itau2,j,k,:) = 1.0D0       
            prim2(irho2,j,k,:) = 1.0D0
            prim2(ivel2_x,j,k,:) = 0.7276D0
            prim2(ivel2_y,j,k,:) = 0.5D0
	      ELSEIF(x2(j) > 0.5D0 .AND. x2(j) < 1.0D0 .AND. y2(k) > 0.5D0 .AND. y2(k) < 1.0D0) THEN
            prim2(itau2,j,k,:) = 0.4D0  
            prim2(irho2,j,k,:) = 0.5313D0 
            prim2(ivel2_x,j,k,:) = 0.0D0
            prim2(ivel2_y,j,k,:) = 0.0D0
	      END IF
	      epsilon2(j,k,:) = prim2(itau2,j,k,:) / prim2(irho2,j,k,:) / (ggas2 - 1.0D0)
      ENDDO
   ENDDO

! 2D Riemann Problem Test with gamma = 1.4
ELSEIF(test_model == 12) THEN

   ggas2 = 1.4D0

   DO k = 1, ny_2
      DO j = 1, nx_2
	      IF(x2(j) > 0.0D0 .AND. x2(j) < 0.5D0 .AND. y2(k) > 0.0D0 .AND. y2(k) < 0.5D0) THEN
            prim2(itau2,j,k,:) = 0.4D0       
            prim2(irho2,j,k,:) = 0.8D0
            prim2(ivel2_x,j,k,:) = 0.1D0
            prim2(ivel2_y,j,k,:) = -0.3D0
	      ELSEIF(x2(j) > 0.5D0 .AND. x2(j) < 1.0D0 .AND. y2(k) > 0.0D0 .AND. y2(k) < 0.5D0) THEN
            prim2(itau2,j,k,:) = 0.4D0       
            prim2(irho2,j,k,:) = 0.5313D0
            prim2(ivel2_x,j,k,:) = 0.1D0
            prim2(ivel2_y,j,k,:) = 0.4276D0
	      ELSEIF(x2(j) > 0.0D0 .AND. x2(j) < 0.5D0 .AND. y2(k) > 0.5D0 .AND. y2(k) < 1.0D0) THEN
            prim2(itau2,j,k,:) = 0.4D0       
            prim2(irho2,j,k,:) = 0.5197D0
            prim2(ivel2_x,j,k,:) = -0.6259D0
            prim2(ivel2_y,j,k,:) = -0.3D0
	      ELSEIF(x2(j) > 0.5D0 .AND. x2(j) < 1.0D0 .AND. y2(k) > 0.5D0 .AND. y2(k) < 1.0D0) THEN
            prim2(itau2,j,k,:) = 1.0D0  
            prim2(irho2,j,k,:) = 1.0D0 
            prim2(ivel2_x,j,k,:) = 0.1D0
            prim2(ivel2_y,j,k,:) = -0.3D0
	      END IF
	      epsilon2(j,k,:) = prim2(itau2,j,k,:) / prim2(irho2,j,k,:) / (ggas2 - 1.0D0)
      ENDDO
   ENDDO

! 2D Riemann Problem Test with gamma = 1.4
ELSEIF(test_model == 13) THEN

   ggas2 = 1.4D0

   DO k = 1, ny_2
      DO j = 1, nx_2
	      IF(x2(j) > 0.0D0 .AND. x2(j) < 0.5D0 .AND. y2(k) > 0.0D0 .AND. y2(k) < 0.5D0) THEN
            prim2(itau2,j,k,:) = 0.4D0       
            prim2(irho2,j,k,:) = 1.0625D0
            prim2(ivel2_x,j,k,:) = 0.0D0
            prim2(ivel2_y,j,k,:) = 0.2145D0
	      ELSEIF(x2(j) > 0.5D0 .AND. x2(j) < 1.0D0 .AND. y2(k) > 0.0D0 .AND. y2(k) < 0.5D0) THEN
            prim2(itau2,j,k,:) = 0.4D0       
            prim2(irho2,j,k,:) = 0.5197D0
            prim2(ivel2_x,j,k,:) = 0.0D0
            prim2(ivel2_y,j,k,:) = -1.1259D0
	      ELSEIF(x2(j) > 0.0D0 .AND. x2(j) < 0.5D0 .AND. y2(k) > 0.5D0 .AND. y2(k) < 1.0D0) THEN
            prim2(itau2,j,k,:) = 1.0D0       
            prim2(irho2,j,k,:) = 2.0D0
            prim2(ivel2_x,j,k,:) = 0.0D0
            prim2(ivel2_y,j,k,:) = -0.3D0
	      ELSEIF(x2(j) > 0.5D0 .AND. x2(j) < 1.0D0 .AND. y2(k) > 0.5D0 .AND. y2(k) < 1.0D0) THEN
            prim2(itau2,j,k,:) = 1.0D0  
            prim2(irho2,j,k,:) = 1.0D0 
            prim2(ivel2_x,j,k,:) = 0.0D0
            prim2(ivel2_y,j,k,:) = -0.4D0
	      END IF
	      epsilon2(j,k,:) = prim2(itau2,j,k,:) / prim2(irho2,j,k,:) / (ggas2 - 1.0D0)
      ENDDO
   ENDDO

ENDIF

! set boundary conditions !
call BOUNDARYP_NM
call BOUNDARY1D_NM (epsilon2,even,part)

! set atmospheric primitive variables !
prim2_a(:) = 0.0D0
eps2_a = 0.0D0
temp2_a = 0.0D0

END SUBROUTINE 