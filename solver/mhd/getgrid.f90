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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Dummy variables
INTEGER :: j, k, l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Grid position of motherboard

! Then, get interface coordinate !
DO j = -3, nx_2 + 3
	xF2(j) = x2_start + DBLE(j)*dx2_ini
END DO
DO j = -3, ny_2 + 3
	yF2(j) = y2_start + DFLOAT(j)*dy2_ini
END DO
DO j = -3, nz_2 + 3
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
! Centroid 

! Along the first direction !
IF(coordinate_flag == 0) THEN
	DO j = -2, nx_2 + 3
		x2cen(j) = x2(j)
	END DO
ELSEIF(coordinate_flag == 1) THEN
	DO j = -2, nx_2 + 3
		x2cen(j) = x2(j) + dx2(j)**2/(12.0d0*x2(j))
	END DO
ELSEIF(coordinate_flag == 2) THEN
	DO j = -2, nx_2 + 3
		x2cen(j) = x2(j) + 2.0d0*x2(j)*dx2(j)**2/(12.0d0*x2(j)**2 + dx2(j)**2)
	END DO
END IF

! Along the second direction !
IF(coordinate_flag == 2) THEN
	DO k = -2, ny_2 + 3
		y2cen(k) = (yF2(k-1)*DCOS(yF2(k-1)) - yF2(k)*DCOS(yF2(k)) + DSIN(yF2(k)) - DSIN(yF2(k-1))) & 
						 / (DCOS(yF2(k-1)) - DCOS(yF2(k)))
	END DO
ELSE
	DO k = -2, ny_2 + 3
		y2cen(k) = y2(k)
	END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! These are dx square and cube

DO j = -1, nx_2 + 3
	dx2_sq(j) = xF2(j)**2 - xF2(j-1)**2
	dx2_cb(j) = xF2(j)**3 - xF2(j-1)**3
	dx2_qd(j) = xF2(j)**4 - xF2(j-1)**4
	x2bar(j) = (2.0D0/3.0D0)*(dx2_cb(j))/(dx2_sq(j))
END DO
DO k = -1, ny_2 + 3
	dcos2(k) = - (DCOS(yF2(k)) - DCOS(yF2(k-1)))
	dsin2(k) = (DSIN(yF2(k)) - DSIN(yF2(k-1)))
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Surface area and differential volume !

! Get the volume, assuming full box simulation ! 
IF(coordinate_flag == 0) THEN
	DO l = -1, nz_2 + 3
		DO k = -1, ny_2 + 3
			DO j = -1, nx_2 + 3
				vol2(j,k,l) = dx2(j)*dy2(k)*dz2(l)
			END DO
		END DO
	END DO
ELSEIF(coordinate_flag == 1) THEN
	DO l = -1, nz_2 + 3
		DO k = -1, ny_2 + 3
			DO j = -1, nx_2 + 3
				vol2(j,k,l) = 0.5D0*(xF2(j)**2 - xF2(j-1)**2)*dy2(k)*dz2(l)
			END DO
		END DO
	END DO
ELSEIF(coordinate_flag == 2) THEN
	DO l = -1, nz_2 + 3
		DO k = -1, ny_2 + 3
			DO j = -1, nx_2 + 3
				vol2(j,k,l) = (xF2(j)**3 - xF2(j-1)**3)*ABS(DCOS(yF2(k)) - DCOS(yF2(k-1)))*dz2(l)/3.0D0
			END DO
		END DO
	END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! volume and surface element cannot be negative

vol2(:,:,:) = ABS(vol2(:,:,:))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE