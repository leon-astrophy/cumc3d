!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine sets up the position of each grid
! Notice if dx is changed, this subroutine needs to be called every time
! This subroutine also finds the local volume to save computing time
! Note that we adopt a coordinate system convention so that 
! Cartesian (x, y, z) Cylindrical (r, z, phi) Spherical (r, theta phi)
!
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE GetGrid
USE definition
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Dummy variables
INTEGER :: j, k, l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Grid position of motherboard

! Then, get interface coordinate !
DO j = -3, nx + 3
	xF(j) = x_start + DBLE(j)*dx_ini
END DO
DO j = -3, ny + 3
	yF(j) = y_start + DBLE(j)*dy_ini
END DO
DO j = -3, nz + 3
	zF(j) = z_start + DBLE(j)*dz_ini
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Or, provide your own grid function if you dislike

CALL CUSTOM_GRID

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! First, assign grid size !
DO j = -2, nx + 3
	dx(j) = xF(j) - xF(j-1)
END DO
DO j = -2, ny + 3
	dy(j) = yF(j) - yF(j-1)
END DO
DO j = -2, nz + 3
	dz(j) = zF(j) - zF(j-1)
END DO

! Now assign distance !
DO j = -2, nx + 3
	x(j) = (xF(j) + xF(j-1))/2
END DO
DO j = -2, ny + 3
	y(j) = (yF(j) + yF(j-1))/2
END DO
DO j = -2, nz + 3
	z(j) = (zF(j) + zF(j-1))/2
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! These are dx square and cube

DO j = -1, nx + 3
	dx_cb(j) = xF(j)**3 - xF(j-1)**3
	xbar(j) = (dx_cb(j)/3.0D0)/(x(j)*dx(j))
END DO
DO k = -1, ny + 3
	dcose(k) = - (DCOS(yF(k)) - DCOS(yF(k-1)))
	dsine(k) = (DSIN(yF(k)) - DSIN(yF(k-1)))
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! sine !

! sine at cell center and surface !
DO k = -2, ny + 3
	sine(k) = DSIN(y(k))
END DO
DO k = -3, ny + 3
	sinf(k) = DSIN(yF(k))
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! fix the angle and sine at coordinate axis !
#ifdef FIXPOLE
IF(coordinate_flag == 2) THEN
	yF(0) = 0.0d0
	yF(ny) = pi
	sinf(0) = 0.0D0
	sinf(ny) = 0.0D0
END IF
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Surface area and differential volume !

! Get the volume, assuming full box simulation ! 
IF(coordinate_flag == 0) THEN
	DO l = -1, nz + 3
		DO k = -1, ny + 3
			DO j = -1, nx + 3
				vol(j,k,l) = dx(j)*dy(k)*dz(l)
			END DO
		END DO
	END DO
ELSEIF(coordinate_flag == 1) THEN
	DO l = -1, nz + 3
		DO k = -1, ny + 3
			DO j = -1, nx + 3
				vol(j,k,l) = 0.5D0*(xF(j)**2 - xF(j-1)**2)*dy(k)*dz(l)
			END DO
		END DO
	END DO
ELSEIF(coordinate_flag == 2) THEN
	DO l = -1, nz + 3
		DO k = -1, ny + 3
			DO j = -1, nx + 3
				vol(j,k,l) = (xF(j)**3 - xF(j-1)**3)*ABS(DCOS(yF(k)) - DCOS(yF(k-1)))*dz(l)/3.0D0
			END DO
		END DO
	END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Special treatment for full box simulation !
IF(NOT(fullx_flag)) THEN
	DO l = -2, nz + 3
		DO k = -2, ny + 3
			DO j = -2, nx + 3
				vol(j,k,l) = vol(j,k,l)*2.0D0
			END DO
		END DO
	END DO
END IF
IF(NOT(fully_flag)) THEN
	DO l = -2, nz + 3
		DO k = -2, ny + 3
			DO j = -2, nx + 3
				vol(j,k,l) = vol(j,k,l)*2.0D0
			END DO
		END DO
	END DO
END IF
IF(NOT(fullz_flag)) THEN
	DO l = -2, nz + 3
		DO k = -2, ny + 3
			DO j = -2, nx + 3
				vol(j,k,l) = vol(j,k,l)*2.0D0
			END DO
		END DO
	END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! volume and surface element cannot be negative

vol(:,:,:) = ABS(vol(:,:,:))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE
