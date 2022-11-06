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
INTEGER :: j, k, l, n_mid

! Real !
REAL (DP) :: dh_x, dh_y, dh_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign !
IF(movinggridnm_flag) THEN
	dh_x = delta2
	dh_y = delta2
	dh_z = delta2
ELSE
	dh_x = dx2
	dh_y = dy2
	dh_z = dz2
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Grid position of motherboard

! Get the x-position
IF(fornax_flag) THEN
    DO j = -2, nx_2 + 3
		xF2(j) = A_fornax * xt * dsinh(dble(j)/xt)
    END DO
    DO j = -2, nx_2 + 3
        x2(j) = 0.5D0*(xF2(j) + xF2(j-1))
        dr2(j) = xF2(j) - xF2(j-1)
    END DO
ELSE
    DO j = -2, nx_2 + 3
        x2(j) = (DBLE(j) - 0.5D0) * dh_x
        xF2(j) = DBLE(j) * dh_x
        dr2(j) = dh_x
    ENDDO
END IF

! Special treatment for full box simulation !
IF(fullx_flag) THEN
	n_mid = nx_2/2
   	DO j = -2, nx_2 + 3
        x2(j) = (DBLE(j - n_mid) - 0.5D0) * dh_x
		xF2(j) = DBLE(j - n_mid)* dh_x
      	dr2(j) = dh_x
   	ENDDO
END IF

! Get the y-position
IF(n_dim > 1) THEN
	DO j = -2, ny_2 + 3
   	 	y2(j) = (DBLE(j) - 0.5D0) * dh_y
    	yF2(j) = DBLE(j)* dh_y
	ENDDO
	IF(fully_flag) THEN
		n_mid = ny_2/2
   		DO j = -2, ny_2 + 3
        	y2(j) = (DBLE(j - n_mid) - 0.5D0) * dh_y
			yF2(j) = DBLE(j - n_mid)* dh_y
   		ENDDO
	END IF
END IF

! Get the z-position
IF(n_dim > 2) THEN
	DO j = -2, nz_2 + 3
   	 	z2(j) = (DBLE(j) - 0.5D0) * dh_z
    	zF2(j) = DBLE(j)* dh_z
	ENDDO
	IF(fullz_flag) THEN
		n_mid = nz_2/2
   		DO j = -2, nz_2 + 3
        	z2(j) = (DBLE(j - n_mid) - 0.5D0) * dh_z
			zF2(j) = DBLE(j - n_mid)* dh_z
   		ENDDO
	END IF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Grid radial distances, and sin, cos ...

! Radial distances !
IF(coordinate_flag == 2) THEN
	DO CONCURRENT (j = -2:nx_2+3, k = -2:ny_2+3, l = -2:nz_2+3)
		rad2(j,k,l) = x2(j)
		IF(n_dim > 1) THEN
 			cos2(j,k,l) = DCOS(y2(k))
			sin2(j,k,l) = DSIN(y2(k))
		END IF
	END DO
ELSEIF(coordinate_flag == 1) THEN
   	DO CONCURRENT (j = -2:nx_2+3, k = -2:ny_2+3, l = -2:nz_2+3)
		rad2(j,k,l) = x2(j)**2 
		IF(n_dim > 1) THEN
			rad2(j,k,l) = rad2(j,k,l) + y2(k)**2
		END IF
		rad2(j,k,l) = SQRT(rad2(j,k,l))
 		cos2(j,k,l) = y2(k)/rad2(j,k,l)
		sin2(j,k,l) = x2(j)/rad2(j,k,l)
	END DO
ELSE
   	DO CONCURRENT (j = -2:nx_2+3, k = -2:ny_2+3, l = -2:nz_2+3)
		rad2(j,k,l) = x2(j)**2
		IF(n_dim > 1) THEN
			rad2(j,k,l) = rad2(j,k,l) + y2(k)**2
		END IF
		IF(n_dim > 2) THEN
			rad2(j,k,l) = rad2(j,k,l) + z2(l)**2
		END IF
		rad2(j,k,l) = SQRT(rad2(j,k,l))
		IF(n_dim > 1) THEN
			sin2(j,k,l) = SQRT(x2(j)**2 + y2(k)**2)/rad2(j,k,l)
		END IF
		IF(n_dim > 2) THEN
			cos2(j,k,l) = z2(l)/rad2(j,k,l)
		END IF
	END DO
END IF

! Get maxmimum distance !
rmax2 = maxval(x2)
IF(n_dim > 1) THEN
	rmax2 = MAX(rmax2, maxval(y2))
END IF
IF(n_dim > 2) THEN
	rmax2 = MAX(rmax2, maxval(z2))
END IF 

! Dimensionless distances !
DO CONCURRENT (j = -2:nx_2+3, k = -2:ny_2+3, l = -2:nz_2+3)
	radbar2(j,k,l) = rad2(j,k,l)/rmax2
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Get the volume, assuming full box simulation ! 
IF(coordinate_flag == 0) THEN
	DO CONCURRENT (j = -2:nx_2+3, k = -2:ny_2+3, l = -2:nz_2+3)
    	vol2(j,k,l) = dh_x * dh_y * dh_z
	END DO
ELSEIF(coordinate_flag == 1) THEN
	DO CONCURRENT (j = -2:nx_2+3, k = -2:ny_2+3, l = -2:nz_2+3)
    	vol2(j,k,l) = 0.5D0 * (xF2 (j)**2 - xF2 (j-1)**2) * dh_y * dh_z
    ENDDO
ELSEIF(coordinate_flag == 2) THEN
    DO CONCURRENT (j = -2:nx_2+3, k = -2:ny_2+3, l = -2:nz_2+3)
        vol2(j,k,l) = (1.0D0/3.0D0) * (xF2 (j)**3 - xF2 (j-1)**3) * ABS(COS(yF2(k)) - COS(yF2(k-1))) * dh_z
    ENDDO
END IF

! Special treatment for full box simulation !
IF(NOT(fullx_flag)) THEN
	DO CONCURRENT (j = -2:nx_2+3, k = -2:ny_2+3, l = -2:nz_2+3)
		vol2(j,k,l) = vol2(j,k,l)*2.0D0
	END DO
END IF
IF(NOT(fully_flag)) THEN
	DO CONCURRENT (j = -2:nx_2+3, k = -2:ny_2+3, l = -2:nz_2+3)
		vol2(j,k,l) = vol2(j,k,l)*2.0D0
	END Do
END IF
IF(NOT(fullz_flag)) THEN
	DO CONCURRENT (j = -2:nx_2+3, k = -2:ny_2+3, l = -2:nz_2+3)
		vol2(j,k,l) = vol2(j,k,l)*2.0D0
	END Do
END IF

! Dimensionless volume !
DO CONCURRENT (j = -2:nx_2+3, k = -2:ny_2+3, l = -2:nz_2+3)
	volbar2(j,k,l) = vol2(j,k,l)/rmax2**3
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do everything the same, but for DM !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE GetGrid_DM
USE definition
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Dummy variables
INTEGER :: j, k, l, n_mid

! Real !
REAL (DP) :: dh_x, dh_y, dh_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign !
IF(movinggriddm_flag) THEN
	dh_x = delta1
	dh_y = delta1
	dh_z = delta1
ELSE
	dh_x = dx1
	dh_y = dy1
	dh_z = dz1
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Grid position of motherboard

! Get the x-position
IF(fornax_flag) THEN
    DO j = -2, nx_1 + 3
		xF1(j) = A_fornax * xt * dsinh(dble(j)/xt)
    END DO
    DO j = -2, nx_1 + 3
        x1(j) = 0.5D0*(xF1(j) + xF1(j-1))
        dr1(j) = xF1(j) - xF1(j-1)
    END DO
ELSE
    DO j = -2, nx_1 + 3
        x1(j) = (DBLE(j) - 0.5D0) * dh_x
        xF1(j) = DBLE(j) * dh_x
        dr1(j) = dh_x
    ENDDO
END IF

! Special treatment for full box simulation !
IF(fullx_flag) THEN
	n_mid = nx_1/2
   	DO j = -2, nx_1 + 3
        x1(j) = (DBLE(j - n_mid) - 0.5D0) * dh_x
		xF1(j) = DBLE(j - n_mid)* dh_x
      	dr1(j) = dh_x
   	ENDDO
END IF

! Get the y-position
IF(n_dim > 1) THEN
	DO j = -2, ny_1 + 3
   	 	y1(j) = (DBLE(j) - 0.5D0) * dh_y
    	yF1(j) = DBLE(j)* dh_y
	ENDDO
	IF(fully_flag) THEN
		n_mid = ny_1/2
   		DO j = -2, ny_1 + 3
        	y1(j) = (DBLE(j - n_mid) - 0.5D0) * dh_y
			yF1(j) = DBLE(j - n_mid)* dh_y
   		ENDDO
	END IF
END IF

! Get the z-position
IF(n_dim > 2) THEN
	DO j = -2, nz_1 + 3
   	 	z1(j) = (DBLE(j) - 0.5D0) * dh_z
    	zF1(j) = DBLE(j)* dh_z
	ENDDO
	IF(fullz_flag) THEN
		n_mid = nz_1/2
   		DO j = -2, nz_1 + 3
        	z1(j) = (DBLE(j - n_mid) - 0.5D0) * dh_z
			zF1(j) = DBLE(j - n_mid)* dh_z
   		ENDDO
	END IF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Grid radial distances, and sin, cos ...

! Radial distances !
IF(coordinate_flag == 2) THEN
	DO CONCURRENT (j = -2:nx_1+3, k = -2:ny_1+3, l = -2:nz_1+3)
		rad1(j,k,l) = x1(j)
		IF(n_dim > 1) THEN
 			cos1(j,k,l) = DCOS(y1(k))
			sin1(j,k,l) = DSIN(y1(k))
		END IF
	END DO
ELSEIF(coordinate_flag == 1) THEN
   	DO CONCURRENT (j = -2:nx_1+3, k = -2:ny_1+3, l = -2:nz_1+3)
		rad1(j,k,l) = x1(j)**2 
		IF(n_dim > 1) THEN
			rad1(j,k,l) = rad1(j,k,l) + y1(k)**2
		END IF
		rad1(j,k,l) = SQRT(rad1(j,k,l))
 		cos1(j,k,l) = y1(k)/rad1(j,k,l)
		sin1(j,k,l) = x1(j)/rad1(j,k,l)
	END DO
ELSE
   	DO CONCURRENT (j = -2:nx_1+3, k = -2:ny_1+3, l = -2:nz_1+3)
		rad1(j,k,l) = x1(j)**2
		IF(n_dim > 1) THEN
			rad1(j,k,l) = rad1(j,k,l) + y1(k)**2
		END IF
		IF(n_dim > 2) THEN
			rad1(j,k,l) = rad1(j,k,l) + z1(l)**2
		END IF
		rad1(j,k,l) = SQRT(rad1(j,k,l))
		IF(n_dim > 1) THEN
			sin1(j,k,l) = SQRT(x1(j)**2 + y1(k)**2)/rad1(j,k,l)
		END IF
		IF(n_dim > 2) THEN
			cos1(j,k,l) = z1(l)/rad1(j,k,l)
		END IF
	END DO
END IF

! Get maxmimum distance !
rmax1 = maxval(x1)
IF(n_dim > 1) THEN
	rmax1 = MAX(rmax1, maxval(y1))
END IF
IF(n_dim > 2) THEN
	rmax1 = MAX(rmax1, maxval(z1))
END IF 

! Dimensionless distances !
DO CONCURRENT (j = -2:nx_1+3, k = -2:ny_1+3, l = -2:nz_1+3)
	radbar1(j,k,l) = rad1(j,k,l)/rmax1
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Get the volume, assuming full box simulation ! 
IF(coordinate_flag == 0) THEN
	DO CONCURRENT (j = -2:nx_1+3, k = -2:ny_1+3, l = -2:nz_1+3)
    	vol1(j,k,l) = dh_x * dh_y * dh_z
	END DO
ELSEIF(coordinate_flag == 1) THEN
	DO CONCURRENT (j = -2:nx_1+3, k = -2:ny_1+3, l = -2:nz_1+3)
    	vol1(j,k,l) = 0.5D0 * (xF1 (j)**2 - xF1 (j-1)**2) * dh_y * dh_z
    ENDDO
ELSEIF(coordinate_flag == 2) THEN
    DO CONCURRENT (j = -2:nx_1+3, k = -2:ny_1+3, l = -2:nz_1+3)
        vol1(j,k,l) = (1.0D0/3.0D0) * (xF1 (j)**3 - xF1 (j-1)**3) * ABS(COS(yF1(k)) - COS(yF1(k-1))) * dh_z
    ENDDO
END IF

! Special treatment for full box simulation !
IF(NOT(fullx_flag)) THEN
	DO CONCURRENT (j = -2:nx_1+3, k = -2:ny_1+3, l = -2:nz_1+3)
		vol1(j,k,l) = vol1(j,k,l)*2.0D0
	END DO
END IF
IF(NOT(fully_flag)) THEN
	DO CONCURRENT (j = -2:nx_1+3, k = -2:ny_1+3, l = -2:nz_1+3)
		vol1(j,k,l) = vol1(j,k,l)*2.0D0
	END Do
END IF
IF(NOT(fullz_flag)) THEN
	DO CONCURRENT (j = -2:nx_1+3, k = -2:ny_1+3, l = -2:nz_1+3)
		vol1(j,k,l) = vol1(j,k,l)*2.0D0
	END Do
END IF

! Dimensionless volume !
DO CONCURRENT (j = -2:nx_1+3, k = -2:ny_1+3, l = -2:nz_1+3)
	volbar1(j,k,l) = vol1(j,k,l)/rmax1**3
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE