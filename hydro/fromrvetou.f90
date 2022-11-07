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
USE DEFINITION
IMPLICIT NONE

! Dummy variables
INTEGER :: i, j, k, l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DM Sector !

! From RVE to U 
IF(DM_flag) then
	DO CONCURRENT (j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1) 
		DO i = imin1, imax1
			cons1(i,j,k,l) = prim1(i,j,k,l)*prim1(irho1,j,k,l)
		END DO
		cons1(irho1,j,k,l) = prim1(irho1,j,k,l)
	END DO
END IF
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NM Sector !

! Convert NM hydro
DO CONCURRENT (j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2) 
	DO i = imin2, imax2
		cons2(i,j,k,l) = prim2(i,j,k,l)*prim2(irho2,j,k,l)
	END DO
	cons2(irho2,j,k,l) = prim2(irho2,j,k,l)
	cons2(itau2,j,k,l) = prim2(irho2,j,k,l)*(epsilon2(j,k,l) + 0.5D0*(prim2(ivel2_x,j,k,l)**2 & 
				  	   		   + prim2(ivel2_y,j,k,l)**2 + prim2(ivel2_z,j,k,l)**2))

	! Dual energy !
	IF(dual_energy) THEN
		cons2(ieps2,j,k,l) = prim2(irho2,j,k,l) * epsilon2(j,k,l)
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
IF(DM_flag) THEN
	DO CONCURRENT (j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1)
		DO i = imin1, imax1
			prim1(i,j,k,l) = cons1(i,j,k,l)/cons1(irho1,j,k,l)
		END DO
		prim1(irho1,j,k,l) = cons1(irho1,j,k,l)
	END DO

	! Copy to boundaries !
   	CALL BOUNDARYP_DM
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For NM sectors !

! Convert the NM hydro
DO CONCURRENT (j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2) 
	DO i = imin2, imax2
		prim2(i,j,k,l) = cons2(i,j,k,l)/cons2(irho2,j,k,l)
	END DO
	prim2(irho2,j,k,l) = cons2(irho2,j,k,l)

	! Determine the epsilon for epsilon equation !
	IF (dual_energy) THEN
		prim2(ieps2,j,k,l) = cons2(ieps2,j,k,l)
		bige_2 = cons2(itau2,j,k,l)/cons2(irho2,j,k,l)
		et_2 = bige_2 - 5.0E-1_DP * (prim2(ivel2_x,j,k,l) ** 2 & 
			 + prim2(ivel2_y,j,k,l) ** 2 + prim2(ivel2_z,j,k,l) ** 2)
		If(et_2 > 1.0D-4*bige_2) THEN
			epsilon2(j,k,l) = et_2 
		ELSE
			epsilon2(j,k,l) = cons2 (ieps2,j,k,l)/cons2(irho2,j,k,l)
		END	IF
	ELSE
		epsilon2(j,k,l) = (cons2(itau2,j,k,l)/cons2(irho2,j,k,l) - 5.0E-1_DP * (prim2(ivel2_x,j,k,l) ** 2 + prim2(ivel2_y,j,k,l) ** 2 + prim2(ivel2_z,j,k,l) ** 2))
	END IF
END DO

! boundary condition !
call BOUNDARYP_NM
call boundary1d_NM (epsilon2,even,part)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE