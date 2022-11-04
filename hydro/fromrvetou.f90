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
INTEGER :: i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DM Sector !

! From RVE to U 
IF(runDM_flag) then
	DO i = imin1, imax1
		cons1(:,:,:,i) = prim1(:,:,:,i)*prim1(:,:,:,irho1)
	END DO
	cons1(:,:,:,irho1) = prim1(:,:,:,irho1)
END IF
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NM Sector !

! Convert NM hydro
DO i = imin2, imax2
	cons2(:,:,:,i) = prim2(:,:,:,i)*prim2(:,:,:,irho2)
END DO
cons2(:,:,:,irho2) = prim2(:,:,:,irho2)
cons2(:,:,:,itau2) = prim2(:,:,:,irho2)*(epsilon2(:,:,:) + 0.5D0*(prim2(:,:,:,ivel2_x)**2 & 
				   + prim2(:,:,:,ivel2_y)**2 + prim2(:,:,:,ivel2_z)**2))

! Dual energy !
IF(dual_energy) THEN
	cons2(:,:,:,ieps2) = prim2(:,:,:,irho2) * epsilon2(:,:,:)
END IF

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
	DO i = imin1, imax1
		prim1(:,:,:,i) = cons1(:,:,:,i)/cons1(:,:,:,irho1)
	END DO
	prim1(:,:,:,irho1) = cons1(:,:,:,irho1)

	! Copy to boundaries !
   	CALL BOUNDARYP_DM
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For NM sectors !

! Convert the NM hydro
DO i = imin2, imax2
	prim2(:,:,:,i) = cons2(:,:,:,i)/cons2(:,:,:,irho2)
END DO
prim2(:,:,:,irho2) = cons2(:,:,:,irho2)

! Determine the epsilon for epsilon equation !
IF (dual_energy) THEN
	prim2(:,:,:,ieps2) = cons2(:,:,:,ieps2)
	DO j = 1, nx_part_2
		DO k = 1, ny_part_2
			DO l = 1, nz_part_2
				et_2 = (cons2(j,k,l,itau2) - 5.0E-1_DP * (prim2(j,k,l,ivel2_x) ** 2 & 
					 + prim2(j,k,l,ivel2_y) ** 2 + prim2(j,k,l,ivel2_z))) / prim2(j,k,l,irho2)
				bige_2 = cons2(j,k,l,itau2) / prim2(j,k,l,irho2)
				If(et_2 > 1.0D-4*bige_2) THEN
					epsilon2(j,k,l) = et_2 
				ELSE
					epsilon2(j,k,l) = cons2 (j,k,l,ieps2) / prim2(j,k,l,irho2)
				END	IF
			END DO
		END DO
	END DO
ELSE
	epsilon2(:,:,:) = (cons2(j,k,l,itau2) - 5.0E-1_DP * (prim2(j,k,l,ivel2_x) ** 2 & 
					+ prim2(j,k,l,ivel2_y) ** 2 + prim2(j,k,l,ivel2_z))) / prim2(j,k,l,irho2)
END IF

! boundary condition !
call BOUNDARYP_NM
call boundary1d_NM (epsilon2,even,part)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE