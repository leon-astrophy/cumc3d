!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine converts the primitive variables to
! conservative variables (or vice versa)
! Written by Leung Shing Chi in 2016
! If you add your own variables in the WENO scheme,
! add your conversion step here.
! This subroutine takes in the U array and conversion mode
! Mode 0: From primitive to conservative
! Mode 1: From conservative to primitive
!
! Here is a reminder in how to add new physics:
! 1. Add your own module that contains the physicsb
! 2. Remind BuildWENO to include your quantity
! 3. Add the conversion here
! 4. Write a section in how to calculate the flux in Spatial
! 5. Add a flag to give signal to the program whenever you use the code
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE FROMRVETOU
USE OMP_LIB
USE DEFINITION
IMPLICIT NONE

! Dummy variables
INTEGER :: i, j, k, l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DM Sector !

! From RVE to U 
IF(runDM_flag) then
	DO CONCURRENT (j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1) 
		DO i = imin1, imax1
			cons1(j,k,l,i) = prim1(j,k,l,i)*prim1(j,k,l,irho1)
		END DO
		cons1(j,k,l,irho1) = prim1(j,k,l,irho1)
	END DO
END IF
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NM Sector !

! Convert NM hydro
DO CONCURRENT (j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2) 
	DO i = imin2, imax2
		cons2(j,k,l,i) = prim2(j,k,l,i)*prim2(j,k,l,irho2)
	END DO
	cons2(j,k,l,irho2) = prim2(j,k,l,irho2)
	cons2(j,k,l,itau2) = prim2(j,k,l,irho2)*(epsilon2(j,k,l) + 0.5D0*(prim2(j,k,l,ivel2_x)**2 & 
				  	   		   + prim2(j,k,l,ivel2_y)**2 + prim2(j,k,l,ivel2_z)**2))

	! Dual energy !
	IF(dual_energy) THEN
		cons2(j,k,l,ieps2) = prim2(j,k,l,irho2) * epsilon2(j,k,l)
	END IF
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Convert back to primitive variables !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FROMUTORVE
USE DEFINITION
IMPLICIT NONE

! Dummy variables
INTEGER :: i, j, k, l

! Dummy variables !
REAL (DP) :: et_2, bige_2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! From U to RVE for DM sectors !

! Convert the DM hydro
IF(runDM_flag) THEN
	DO CONCURRENT (j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1)
		DO i = imin1, imax1
			prim1(j,k,l,i) = cons1(j,k,l,i)/cons1(j,k,l,irho1)
		END DO
		prim1(j,k,l,irho1) = cons1(j,k,l,irho1)
	END DO

	! Copy to boundaries !
   	CALL BOUNDARYP_DM
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For NM sectors !

! Convert the NM hydro
DO CONCURRENT (j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2) 
	DO i = imin2, imax2
		prim2(j,k,l,i) = cons2(j,k,l,i)/cons2(j,k,l,irho2)
	END DO
	prim2(j,k,l,irho2) = cons2(j,k,l,irho2)

	! Determine the epsilon for epsilon equation !
	IF (dual_energy) THEN
		prim2(j,k,l,ieps2) = cons2(j,k,l,ieps2)
		bige_2 = cons2(j,k,l,itau2)/cons2(j,k,l,irho2)
		et_2 = bige_2 - 5.0E-1_DP * (prim2(j,k,l,ivel2_x) ** 2 & 
			 + prim2(j,k,l,ivel2_y) ** 2 + prim2(j,k,l,ivel2_z) ** 2)
		If(et_2 > 1.0D-4*bige_2) THEN
			epsilon2(j,k,l) = et_2 
		ELSE
			epsilon2(j,k,l) = cons2 (j,k,l,ieps2)/cons2(j,k,l,irho2)
		END	IF
	ELSE
		epsilon2(j,k,l) = (cons2(j,k,l,itau2)/cons2(j,k,l,irho2) - 5.0E-1_DP * (prim2(j,k,l,ivel2_x) ** 2 + prim2(j,k,l,ivel2_y) ** 2 + prim2(j,k,l,ivel2_z) ** 2))
	END IF
END DO

! boundary condition !
call BOUNDARYP_NM
call boundary1d_NM (epsilon2,even,part)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE