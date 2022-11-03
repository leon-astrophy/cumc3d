!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This files contain all the riemann solvers available for !
! simulating hydrodynamics				   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE RIEMANN_MODULE
USE DEFINITION
IMPLICIT NONE

! the alpha in the LF flux !
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: alpha1_r, alpha1_z
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: alpha2_r, alpha2_z

! Left and right hydro-states for DM/NM !
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: eps1R, eps1L 
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: eps2R, eps2L

! DM moving grid !
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: vf1rR, vf1rL
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: vf1zR, vf1zL

! NM moving grid !
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: vf2rR, vf2rL
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: vf2zR, vf2zL

! Pressure !
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: p1L, p1R
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: p2L, p2R

! Speed of sound !
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: cs1L, cs1R
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: cs2L, cs2R

! Lapse !
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: ap1L, ap1R
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: ap2L, ap2R

! Left and right fluxes, conserved quantity !
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:) :: prim1, primL1, primR1
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:) :: prim2, primL2, primR2
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:) :: fluxL1, fluxR1, uL1, uR1
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:) :: fluxL2, fluxR2, uL2, uR2

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Assign left and right states and fluxes for riemann problem !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BUILDRIEMANN
USE DEFINITION 
IMPLICIT NONE

! Left and right fluxes, conserved quantity !
IF(RUNDM_flag == 1) THEN
	ALLOCATE (alpha1_r(-4 : length_step_z_1 + 5, imin1:imax1))
	ALLOCATE (alpha1_z(-4 : length_step_r_1 + 5, imin1:imax1))
	ALLOCATE(p1L(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE(p1R(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE(eps1L(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE(eps1R(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE(cs1L(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE(cs1R(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE(fluxL1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5, imin1 : imax1))
	ALLOCATE(fluxR1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5, imin1 : imax1))
	ALLOCATE(uL1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5, imin1 : imax1))
	ALLOCATE(uR1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5, imin1 : imax1))
	ALLOCATE(prim1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5, imin1 : imax1))
	ALLOCATE(primL1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5, imin1 : imax1))
	ALLOCATE(primR1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5, imin1 : imax1))
	IF(lapse_flag == 1) THEN
		ALLOCATE (ap1L(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
		ALLOCATE (ap1R(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	END IF
	IF(movinggriddm_flag == 1 .OR. found_movinggriddm_flag == 1) THEN
		ALLOCATE (vf1rL(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
		ALLOCATE (vf1rR(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
		ALLOCATE (vf1zL(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
		ALLOCATE (vf1zR(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	END IF
END IF

! NM !
ALLOCATE(alpha2_r(-4 : length_step_z_2 + 5, imin2:imax2))
ALLOCATE(alpha2_z(-4 : length_step_r_2 + 5, imin2:imax2))
ALLOCATE(p2L(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE(p2R(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE(eps2L(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE(eps2R(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE(cs2L(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE(cs2R(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE(fluxL2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5, imin2 : imax2))
ALLOCATE(fluxR2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5, imin2 : imax2))
ALLOCATE(uL2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5, imin2 : imax2))
ALLOCATE(uR2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5, imin2 : imax2))
ALLOCATE(prim2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5, imin2 : imax2))
ALLOCATE(primL2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5, imin2 : imax2))
ALLOCATE(primR2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5, imin2 : imax2))
IF(lapse_flag == 1) THEN
	ALLOCATE (ap2L(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
	ALLOCATE (ap2R(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
END IF

! Moving grid !
IF(movinggridnm_flag == 1 .OR. found_movinggridnm_flag == 1) THEN
	ALLOCATE (vf2rL(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
	ALLOCATE (vf2rR(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
	ALLOCATE (vf2zL(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
	ALLOCATE (vf2zR(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Assign left and right states and fluxes for riemann problem !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DESTROYRIEMANN
USE DEFINITION 
IMPLICIT NONE

! Left and right fluxes, conserved quantity !
IF(RUNDM_flag == 1) THEN
	DEALLOCATE (alpha1_r)
	DEALLOCATE (alpha1_z)
	DEALLOCATE(p1L)
	DEALLOCATE(p1R)
	DEALLOCATE(eps1L)
	DEALLOCATE(eps1R)
	DEALLOCATE(cs1L)
	DEALLOCATE(cs1R)
	DEALLOCATE(fluxL1)
	DEALLOCATE(fluxR1)
	DEALLOCATE(uL1)
	DEALLOCATE(uR1)
	DEALLOCATE(prim1)
	DEALLOCATE(primL1)
	DEALLOCATE(primR1)
	IF(lapse_flag == 1) THEN
		DEALLOCATE (ap1L)
		DEALLOCATE (ap1R)
	END IF
	IF(movinggriddm_flag == 1 .OR. found_movinggriddm_flag == 1) THEN
		DEALLOCATE (vf1rL)
		DEALLOCATE (vf1rR)
		DEALLOCATE (vf1zL)
		DEALLOCATE (vf1zR)
	END IF
END IF

! NM !
DEALLOCATE(alpha2_r)
DEALLOCATE(alpha2_z)
DEALLOCATE(p2L)
DEALLOCATE(p2R)
DEALLOCATE(eps2L)
DEALLOCATE(eps2R)
DEALLOCATE(cs2L)
DEALLOCATE(cs2R)
DEALLOCATE(fluxL2)
DEALLOCATE(fluxR2)
DEALLOCATE(uL2)
DEALLOCATE(uR2)
DEALLOCATE(prim2)
DEALLOCATE(primL2)
DEALLOCATE(primR2)
IF(lapse_flag == 1) THEN
	DEALLOCATE (ap2L)
	DEALLOCATE (ap2R)
END IF

! Moving grid !
IF(movinggridnm_flag == 1 .OR. found_movinggridnm_flag == 1) THEN
	DEALLOCATE (vf2rL)
	DEALLOCATE (vf2rR)
	DEALLOCATE (vf2zL)
	DEALLOCATE (vf2zR)
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The LF flux, see for example. shu et al. for the original WENO paper !					    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LFDM_R (flux_out)
USE DEFINITION 
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k

! Output !
REAL (DP), INTENT (OUT), DIMENSION (-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5, imin1 : imax1) :: flux_out

! Temporal arrays !
REAL (DP), DIMENSION (-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5) :: temporal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For DM !
DO j = 0, length_step_r_part_1
	flux_out(j,:,:) = 0.5D0 * (fluxL1 (j,:,:) + fluxR1 (j,:,:) - alpha1_r(:,:) * (uR1 (j,:,:) - uL1 (j,:,:)))
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The LF flux, see for example. shu et al. for the original WENO paper !					    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LFDM_Z (flux_out)
USE DEFINITION 
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k

! output !
REAL (DP), INTENT (OUT), DIMENSION (-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5, imin1 : imax1) :: flux_out

! Temporal arrays !
REAL (DP), DIMENSION (-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5) :: temporal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For DM !
DO k = length_step_z_min_part_1 - 1, length_step_z_part_1
	flux_out(:,k,:) = 0.5D0 * (fluxL1 (:,k,:) + fluxR1 (:,k,:) - alpha1_z(:,:) * (uR1 (:,k,:) - uL1 (:,k,:)))
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The LF flux, see for example. shu et al. for the original WENO paper !					    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LFNM_R (flux_out)
USE DEFINITION 
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k

! Output !
REAL (DP), INTENT (OUT), DIMENSION (-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5, imin2 : imax2) :: flux_out

! Temporal arrays !
REAL (DP), DIMENSION (-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5) :: temporal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For NM !
DO j = 0, length_step_r_part_2
	flux_out(j,:,:) = 0.5D0 * (fluxL2 (j,:,:) + fluxR2 (j,:,:) - alpha2_r(:,:) * (uR2 (j,:,:) - uL2 (j,:,:)))
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The LF flux, see for example. shu et al. for the original WENO paper !					    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LFNM_Z (flux_out)
USE DEFINITION 
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k

! output !
REAL (DP), INTENT (OUT), DIMENSION (-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5, imin2 : imax2) :: flux_out

! Temporal arrays !
REAL (DP), DIMENSION (-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5) :: temporal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For NM !
DO k = length_step_z_min_part_2 - 1, length_step_z_part_2
	flux_out(:,k,:) = 0.5D0 * (fluxL2 (:,k,:) + fluxR2 (:,k,:) - alpha2_z(:,:) * (uR2 (:,k,:) - uL2 (:,k,:)))
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The HLL riemann solver, see Toro. 2009 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE HLLNM_R (flux_out)
USE DEFINITION 
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! output !
REAL (DP), INTENT (OUT), DIMENSION (-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5, imin2 : imax2) :: flux_out

! some local array !
REAL (DP) :: sL, sR

! For more general signal speed !
REAL (DP) :: ubar

! For more general signal speed !
REAL (DP) :: dbar, dbar2, neta2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Compute fluxes !
DO j = 0, length_step_r_part_2 + 1
	DO k = length_step_z_min_part_2 - 1, length_step_z_part_2 + 1
		ubar = compute_roe(primL2(j,k,ivel2_r), primR2(j,k,ivel2_r), primL2(j,k,irho2), primR2(j,k,irho2))
		neta2 = 0.5D0*sqrt(primL2(j,k,irho2)*primR2(j,k,irho2))/(sqrt(primL2(j,k,irho2)) + sqrt(primR2(j,k,irho2)))**2
		dbar2 = compute_roe(cs2L(j,k)**2,cs2R(j,k)**2,primL2(j,k,irho2),primR2(j,k,irho2)) + neta2*(primR2(j,k,ivel2_r) - primL2(j,k,ivel2_r))**2
		dbar = sqrt(dbar2)
		sL = min(primL2(j,k,ivel2_r) - cs2L(j,k), ubar - dbar)
		sR = max(primR2(j,k,ivel2_r) + cs2R(j,k), ubar + dbar)
		IF(sL >= 0.0D0) THEN
			flux_out(j,k,:) = fluxL2(j,k,:)
		ELSEIF(sL < 0.0D0 .AND. sR > 0.0D0) THEN
			DO i = imin2, imax2
				flux_out(j,k,i) = compute_fluxhll(fluxL2(j,k,i),fluxR2(j,k,i),uL2(j,k,i),uR2(j,k,i),sL,sR)
			END DO
		ELSEIF(sR <= 0.0D0) THEN
			flux_out(j,k,:) = fluxR2(j,k,:)
		END IF
	END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

	real(DP) function compute_fluxhll(yl,yr,xl,xr,sl,sr)
	implicit none
	real(DP) :: yl, yr, xl, xr, sl, sr
	compute_fluxhll = (sr*yl-sl*yr+sl*sr*(xr-xl))/(sr-sl)
	end function

	real(DP) function compute_roe(xl,xr,rhol,rhor)
	implicit none
	real(DP) :: xl,xr,rhol,rhor
	compute_roe = (sqrt(rhol)*xl + sqrt(rhor)*xr)/(sqrt(rhol) + sqrt(rhor))
	end function

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The HLL riemann solver, see Toro. 2009 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE HLLNM_Z (flux_out)
USE DEFINITION 
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! output !
REAL (DP), INTENT (OUT), DIMENSION (-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5, imin2 : imax2) :: flux_out

! some local array !
REAL (DP) :: sL, sR

! For more general signal speed !
REAL (DP) :: ubar

! For more general signal speed !
REAL (DP) :: dbar, dbar2, neta2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Compute fluxes !
DO j = 0, length_step_r_part_2 + 1
	DO k = length_step_z_min_part_2 - 1, length_step_z_part_2 + 1
		ubar = compute_roe(primL2(j,k,ivel2_z), primR2(j,k,ivel2_z), primL2(j,k,irho2), primR2(j,k,irho2))
		neta2 = 0.5D0*sqrt(primL2(j,k,irho2)*primR2(j,k,irho2))/(sqrt(primL2(j,k,irho2)) + sqrt(primR2(j,k,irho2)))**2
		dbar2 = compute_roe(cs2L(j,k)**2,cs2R(j,k)**2,primL2(j,k,irho2),primR2(j,k,irho2)) + neta2*(primR2(j,k,ivel2_z) - primL2(j,k,ivel2_z))**2
		dbar = sqrt(dbar2)
		sL = min(primL2(j,k,ivel2_z) - cs2L(j,k), ubar - dbar)
		sR = max(primR2(j,k,ivel2_z) + cs2R(j,k), ubar + dbar)
		IF(sL >= 0.0D0) THEN
			flux_out(j,k,:) = fluxL2(j,k,:)
		ELSEIF(sL < 0.0D0 .AND. sR > 0.0D0) THEN
			DO i = imin2, imax2
				flux_out(j,k,i) = compute_fluxhll(fluxL2(j,k,i),fluxR2(j,k,i),uL2(j,k,i),uR2(j,k,i),sL,sR)
			END DO
		ELSEIF(sR <= 0.0D0) THEN
			flux_out(j,k,:) = fluxR2(j,k,:)
		END IF
	END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

	real(DP) function compute_fluxhll(yl,yr,xl,xr,sl,sr)
	implicit none
	real(DP) :: yl, yr, xl, xr, sl, sr
	compute_fluxhll = (sr*yl-sl*yr+sl*sr*(xr-xl))/(sr-sl)
	end function

	real(DP) function compute_roe(xl,xr,rhol,rhor)
	implicit none
	real(DP) :: xl,xr,rhol,rhor
	compute_roe = (sqrt(rhol)*xl + sqrt(rhor)*xr)/(sqrt(rhol) + sqrt(rhor))
	end function

END SUBROUTINE

END MODULE