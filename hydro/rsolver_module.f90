!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This files contain all the riemann solvers available for !
! simulating hydrodynamics				   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE RIEMANN_MODULE
USE DEFINITION
IMPLICIT NONE

! the alpha in the LF flux !
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:) :: alpha1_x, alpha1_y, alpha1_z
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:) :: alpha2_x, alpha2_y, alpha2_z

! Left and right hydro-states for DM/NM !
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:) :: eps1R, eps1L 
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:) :: eps2R, eps2L

! DM moving grid !
REAL (DP), ALLOCATABLE, DIMENSION(:) :: vf1xR, vf1xL
REAL (DP), ALLOCATABLE, DIMENSION(:) :: vf1yR, vf1yL
REAL (DP), ALLOCATABLE, DIMENSION(:) :: vf1zR, vf1zL

! NM moving grid !
REAL (DP), ALLOCATABLE, DIMENSION(:) :: vf2xR, vf2xL
REAL (DP), ALLOCATABLE, DIMENSION(:) :: vf2yR, vf2yL
REAL (DP), ALLOCATABLE, DIMENSION(:) :: vf2zR, vf2zL

! Pressure !
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:) :: p1L, p1R
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:) :: p2L, p2R

! Speed of sound !
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:) :: cs1L, cs1R
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:) :: cs2L, cs2R

! Lapse !
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:) :: ap1L, ap1R
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:) :: ap2L, ap2R

! Left and right fluxes, conserved quantity !
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:,:) :: primL1, primR1
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:,:) :: primL2, primR2
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:,:) :: fluxL1, fluxR1, uL1, uR1
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:,:) :: fluxL2, fluxR2, uL2, uR2

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Assign left and right states and fluxes for riemann problem !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BUILDRIEMANN
USE DEFINITION 
IMPLICIT NONE

! Left and right fluxes, conserved quantity !
IF(RUNDM_flag) THEN
	ALLOCATE(alpha1_x(-2:ny_1+3,-2:nz_1+3,imin1:imax1))
	ALLOCATE(alpha1_y(-2:nx_1+3,-2:nz_1+3,imin1:imax1))
	ALLOCATE(alpha1_z(-2:nx_1+3,-2:ny_1+3,imin1:imax1))

	ALLOCATE(p1L(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
	ALLOCATE(p1R(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
	ALLOCATE(eps1L(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
	ALLOCATE(eps1R(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
	ALLOCATE(cs1L(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
	ALLOCATE(cs1R(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
	ALLOCATE(fluxL1(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3,imin1:imax1))
	ALLOCATE(fluxR1(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3,imin1:imax1))
	ALLOCATE(uL1(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3,imin1:imax1))
	ALLOCATE(uR1(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3,imin1:imax1))
	ALLOCATE(primL1(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3,imin1:imax1))
	ALLOCATE(primR1(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3,imin1:imax1))
	
	IF(lapse_flag) THEN
		ALLOCATE (ap1L(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3,imin1:imax1))
		ALLOCATE (ap1R(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3,imin1:imax1))
	END IF

	IF(movinggriddm_flag) THEN
		ALLOCATE (vf1xL(-2:nx_1+3))
		ALLOCATE (vf1xR(-2:nx_1+3))
		ALLOCATE (vf1yL(-2:ny_1+3))
		ALLOCATE (vf1yR(-2:ny_1+3))
		ALLOCATE (vf1zL(-2:nz_1+3))
		ALLOCATE (vf1zR(-2:nz_1+3))
	END IF
END IF

! NM !
ALLOCATE(alpha2_x(-2:ny_2+3,-2:nz_2+3,imin2:imax2))
ALLOCATE(alpha2_y(-2:nx_2+3,-2:nz_2+3,imin2:imax2))
ALLOCATE(alpha2_z(-2:nx_2+3,-2:ny_2+3,imin2:imax2))

ALLOCATE(p2L(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(p2R(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(eps2L(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(eps2R(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(cs2L(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(cs2R(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(fluxL2(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3,imin2:imax2))
ALLOCATE(fluxR2(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3,imin2:imax2))
ALLOCATE(uL2(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3,imin2:imax2))
ALLOCATE(uR2(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3,imin2:imax2))
ALLOCATE(primL2(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3,imin2:imax2))
ALLOCATE(primR2(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3,imin2:imax2))

IF(lapse_flag) THEN
	ALLOCATE (ap2L(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
	ALLOCATE (ap2R(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
END IF

! Moving grid !
IF(movinggridnm_flag) THEN
	ALLOCATE (vf2xL(-2:nx_2+3))
	ALLOCATE (vf2xR(-2:nx_2+3))
	ALLOCATE (vf2yL(-2:ny_2+3))
	ALLOCATE (vf2yR(-2:ny_2+3))
	ALLOCATE (vf2zL(-2:nz_2+3))
	ALLOCATE (vf2zR(-2:nz_2+3))
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The LF flux, see for example. shu et al. for the original WENO paper !					    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LFDM_X (flux_out)
USE DEFINITION 
IMPLICIT NONE

! Integer !
INTEGER :: j

! Output !
REAL (DP), INTENT (OUT), DIMENSION (-2:nx_1+3,-2:ny_1+3,-2:nz_1+3,imin1:imax1) :: flux_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For DM !
DO j = nx_min_1 - 1, nx_part_1
	flux_out(j,:,:,:) = 0.5D0 * (fluxL1 (j,:,:.:) + fluxR1 (j,:,:.:) - alpha1_x(:,:,:) * (uR1 (j,:,:,:) - uL1 (j,:,:,:)))
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The LF flux, see for example. shu et al. for the original WENO paper !					    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LFDM_Y (flux_out)
USE DEFINITION 
IMPLICIT NONE

! Integer !
INTEGER :: j

! Output !
REAL (DP), INTENT (OUT), DIMENSION (-2:nx_1+3,-2:ny_1+3,-2:nz_1+3,imin1:imax1) :: flux_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For DM !
DO j = ny_min_1 - 1, ny_part_1
	flux_out(:,j,:,:) = 0.5D0 * (fluxL1 (:,j,:.:) + fluxR1 (:,j,:.:) - alpha1_y(:,:,:) * (uR1 (:,j,:,:) - uL1 (:,j,:,:)))
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
INTEGER :: j

! Output !
REAL (DP), INTENT (OUT), DIMENSION (-2:nx_1+3,-2:ny_1+3,-2:nz_1+3,imin1:imax1) :: flux_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For DM !
DO j = nz_min_1 - 1, nz_part_1
	flux_out(:,:,j,:) = 0.5D0 * (fluxL1 (:,:,j,:) + fluxR1 (:,:,j,:) - alpha1_z(:,:,:) * (uR1 (:,:,j,:) - uL1 (:,:,j,:)))
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The LF flux, see for example. shu et al. for the original WENO paper !					    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LFNM_X (flux_out)
USE DEFINITION 
IMPLICIT NONE

! Integer !
INTEGER :: j

! Output !
REAL (DP), INTENT (OUT), DIMENSION (-2:nx_2+3,-2:ny_2+3,-2:nz_2+3,imin2:imax2) :: flux_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For DM !
DO j = nx_min_2 - 1, nx_part_2
	flux_out(j,:,:,:) = 0.5D0 * (fluxL2 (j,:,:.:) + fluxR2 (j,:,:.:) - alpha2_x(:,:,:) * (uR2 (j,:,:,:) - uL2 (j,:,:,:)))
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The LF flux, see for example. shu et al. for the original WENO paper !					    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LFNM_Y (flux_out)
USE DEFINITION 
IMPLICIT NONE

! Integer !
INTEGER :: j

! Output !
REAL (DP), INTENT (OUT), DIMENSION (-2:nx_2+3,-2:ny_2+3,-2:nz_2+3,imin2:imax2) :: flux_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For DM !
DO j = ny_min_2 - 1, ny_part_2
	flux_out(:,j,:,:) = 0.5D0 * (fluxL2 (:,j,:,:) + fluxR2 (:,j,:,:) - alpha2_y(:,:,:) * (uR2 (:,j,:,:) - uL2 (:,j,:,:)))
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
INTEGER :: j

! Output !
REAL (DP), INTENT (OUT), DIMENSION (-2:nx_2+3,-2:ny_2+3,-2:nz_2+3,imin2:imax2) :: flux_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For DM !
DO j = nz_min_2 - 1, nz_part_2
	flux_out(:,:,j,:) = 0.5D0 * (fluxL2 (:,:,j,:) + fluxR2 (:,:,j,:) - alpha2_z(:,:,j,:) * (uR2 (:,:,j,:) - uL2 (:,:,j,:)))
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

