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

SUBROUTINE GetGrid_NM
USE definition
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Dummy variables
INTEGER :: j, k, l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Grid position of motherboard

! Then, get interface coordinate !
DO j = -2, nx_2 + 3
	xF2(j) = x2_start + DBLE(j)*dx2_ini
END DO
DO j = -2, ny_2 + 3
	yF2(j) = y2_start + DBLE(j)*dy2_ini
END DO
DO j = -2, nz_2 + 3
	zF2(j) = z2_start + DBLE(j)*dz2_ini
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Or, provide your own grid function if you dislike

CALL CUSTOM_GRID

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! First, assign grid size !
DO j = -2, nx_2 + 3
	dx2(j) = xF2(j) - xF2(j-1)
END DO
DO j = -2, ny_2 + 3
	dy2(j) = yF2(j) - yF2(j-1)
END DO
DO j = -2, nz_2 + 3
	dz2(j) = zF2(j) - zF2(j-1)
END DO

! Now assign distance !
DO j = -2, nx_2 + 3
	x2(j) = (xF2(j) + xF2(j-1))/2
END DO
DO j = -2, ny_2 + 3
	y2(j) = (yF2(j) + yF2(j-1))/2
END DO
DO j = -2, nz_2 + 3
	z2(j) = (zF2(j) + zF2(j-1))/2
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Geometric factor for cylindrical/spherical coordinates

! Initialize !
geom2_x = 0.0D0
geom2_y = 0.0D0

IF(coordinate_flag == 1) THEN
	geom2_x(ivel2_y) = 1.0D0
ELSEIF(coordinate_flag == 2) THEN
	geom2_x(ivel2_y) = 1.0D0
	geom2_x(ivel2_z) = 1.0D0
	geom2_y(ivel2_z) = 1.0D0
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Surface area and differential volume !

! Get the volume, assuming full box simulation ! 
IF(coordinate_flag == 0) THEN
	DO l = -2, nz_2 + 3
		DO k = -2, ny_2 + 3
			DO j = -2, nx_2 + 3
				vol2(j,k,l) = dx2(j)*dy2(k)*dz2(l)
			END DO
		END DO
	END DO
ELSEIF(coordinate_flag == 1) THEN
	DO l = -2, nz_2 + 3
		DO k = -2, ny_2 + 3
			DO j = -2, nx_2 + 3
				vol2(j,k,l) = x2(j)*dx2(j)*dy2(k)*dz2(l)
			END DO
		END DO
	END DO
ELSEIF(coordinate_flag == 2) THEN
	DO l = -2, nz_2 + 3
		DO k = -2, ny_2 + 3
			DO j = -2, nx_2 + 3
				vol2(j,k,l) = x2(j)**2*SIN(y2(k))*dx2(j)*dy2(k)*dz2(l)
			END DO
		END DO
	END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Special treatment for full box simulation !
IF(NOT(fullx_flag)) THEN
	DO l = -2, nz_2 + 3
		DO k = -2, ny_2 + 3
			DO j = -2, nx_2 + 3
				vol2(j,k,l) = vol2(j,k,l)*2.0D0
			END DO
		END DO
	END DO
END IF
IF(NOT(fully_flag)) THEN
	DO l = -2, nz_2 + 3
		DO k = -2, ny_2 + 3
			DO j = -2, nx_2 + 3
				vol2(j,k,l) = vol2(j,k,l)*2.0D0
			END DO
		END DO
	END DO
END IF
IF(NOT(fullz_flag)) THEN
	DO l = -2, nz_2 + 3
		DO k = -2, ny_2 + 3
			DO j = -2, nx_2 + 3
				vol2(j,k,l) = vol2(j,k,l)*2.0D0
			END DO
		END DO
	END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! volume and line element cannot be negative

vol2(:,:,:) = ABS(vol2(:,:,:))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE