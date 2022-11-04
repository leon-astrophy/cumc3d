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

! Get the y-position
DO j = -2, ny_2 + 3
    y2(j) = (DBLE(j) - 0.5D0) * dh_y
    yF2(j) = DBLE(j)* dh_y
ENDDO

! Get the z-position
DO j = -2, nz_2 + 3
    z2(j) = (DBLE(j) - 0.5D0) * dh_z
    zF2(j) = DBLE(j)* dh_z
ENDDO

! Special treatment for full box simulation !
IF(fullx_flag) THEN
	n_mid = nx_2/2
   	DO j = -2, nx_2 + 3
        x2(j) = (DBLE(j - n_mid) - 0.5D0) * dh_x
		xF2(j) = DBLE(j - n_mid)* dh_x
      	dr2(j) = dh_x
   	ENDDO
END IF
IF(fully_flag) THEN
	n_mid = ny_2/2
   	DO j = -2, ny_2 + 3
        y2(j) = (DBLE(j - n_mid) - 0.5D0) * dh_y
		yF2(j) = DBLE(j - n_mid)* dh_y
   	ENDDO
END IF
IF(fullz_flag) THEN
	n_mid = nz_2/2
   	DO j = -2, nz_2 + 3
        z2(j) = (DBLE(j - n_mid) - 0.5D0) * dh_z
		zF2(j) = DBLE(j - n_mid)* dh_z
   	ENDDO
END IF

! Radial distances !
IF(coordinate_flag == 2) THEN
   	DO j = -2, nx_2 + 3
		rad2(j,:,:) = x2(j)
	END DO
	DO k = -2, ny_2 + 3
 		cos2(:,k,:) = DCOS(y2(k))
		sin2(:,k,:) = DSIN(y2(k))
   	END DO
ELSEIF(coordinate_flag == 1) THEN
   	DO j = -2, nx_2 + 3
		DO k = -2, ny_2 + 3
			rad2(j,k,:) = SQRT(x2(j)**2 + y2(k)**2)
 			cos2(j,k,:) = y2(k)/rad2(j,k,:)
			sin2(j,k,:) = x2(k)/rad2(j,k,:)
		END DO
   	END DO
ELSE
   	DO j = -2, nx_2 + 3
		DO k = -2, ny_2 + 3
			DO l = -2, nz_2 + 3
				rad2(j,k,l) = SQRT(x2(j)**2 + y2(k)**2 + z2(l)**2)
 				cos2(j,k,l) = z2(k)/rad2(j,k,l)
				sin2(j,k,l) = SQRT(x2(j)**2 + y2(k)**2)/rad2(j,k,l)
			END DO
		END DO
   	END DO
END IF

! Get maxmimum distance !
rmax2 = max(maxval(x2), maxval(y2) , maxval(z2))

! Dimensionless distances !
radbar2(:,:,:) = rad2(:,:,:)/rmax2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Get the volume, assuming full box simulation ! 
IF(coordinate_flag == 0) THEN
    vol2(:,:,:) = dh_x * dh_y * dh_z
ELSEIF(coordinate_flag == 1) THEN
   	DO j = -2, nx_2 + 3
    	vol2(j,:,:) = 0.5D0 * (xF2 (j)**2 - xF2 (j-1)**2) * dh_y * dh_z
    ENDDO
ELSEIF(coordinate_flag == 2) THEN
    DO j = -2, nx_2 + 3
	 	DO k = -2, ny_2 + 3
           vol2(j,k,:) = (1.0D0/3.0D0) * (xF2 (j)**3 - xF2 (j-1)**3) * ABS(COS(yF2(k)) - COS(yF2(k-1))) * dh_z
	 	END DO
    ENDDO
END IF

! Special treatment for full box simulation !
IF(NOT(fullx_flag)) THEN
	vol2(:,:,:) = vol2(:,:,:)*2.0D0
END IF
IF(NOT(fully_flag)) THEN
	vol2(:,:,:) = vol2(:,:,:)*2.0D0
END IF
IF(NOT(fullz_flag)) THEN
	vol2(:,:,:) = vol2(:,:,:)*2.0D0
END IF

! Dimensionless volume !
volbar2(:,:,:) = vol2(:,:,:)/rmax2**3

! Set minimum and maximum domain !
nx_min_2 = 1
nx_part_2 = nx_2
ny_min_2 = 1
ny_part_2 = ny_2
nz_min_2 = 1
nz_part_2 = nz_2

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
IF(movinggridnm_flag) THEN
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

! Get the y-position
DO j = -2, ny_1 + 3
    y1(j) = (DBLE(j) - 0.5D0) * dh_y
    yF1(j) = DBLE(j)* dh_y
ENDDO

! Get the z-position
DO j = -2, nz_1 + 3
    z1(j) = (DBLE(j) - 0.5D0) * dh_z
    zF1(j) = DBLE(j)* dh_z
ENDDO

! Special treatment for full box simulation !
IF(fullx_flag) THEN
	n_mid = nx_1/2
   	DO j = -2, nx_1 + 3
        x1(j) = (DBLE(j - n_mid) - 0.5D0) * dh_x
		xF1(j) = DBLE(j - n_mid)* dh_x
      	dr1(j) = dh_x
   	ENDDO
END IF
IF(fully_flag) THEN
	n_mid = ny_1/2
   	DO j = -2, ny_1 + 3
        y1(j) = (DBLE(j - n_mid) - 0.5D0) * dh_y
		yF1(j) = DBLE(j - n_mid)* dh_y
   	ENDDO
END IF
IF(fullz_flag) THEN
	n_mid = nz_1/2
   	DO j = -2, nz_1 + 3
        z1(j) = (DBLE(j - n_mid) - 0.5D0) * dh_z
		zF1(j) = DBLE(j - n_mid)* dh_z
   	ENDDO
END IF

! Radial distances !
IF(coordinate_flag == 2) THEN
   	DO j = -2, nx_1 + 3
		rad1(j,:,:) = x1(j)
	END DO
	DO k = -2, ny_1 + 3
 		cos1(:,k,:) = DCOS(y1(k))
		sin1(:,k,:) = DSIN(y1(k))
   	END DO
ELSEIF(coordinate_flag == 1) THEN
   	DO j = -2, nx_1 + 3
		DO k = -2, ny_1 + 3
			rad1(j,k,:) = SQRT(x1(j)**2 + y1(k)**2)
 			cos1(j,k,:) = y1(k)/rad1(j,k,:)
			sin1(j,k,:) = x1(k)/rad1(j,k,:)
		END DO
   	END DO
ELSE
   	DO j = -2, nx_1 + 3
		DO k = -2, ny_1 + 3
			DO l = -2, nz_1 + 3
				rad1(j,k,l) = SQRT(x1(j)**2 + y1(k)**2 + z1(l)**2)
 				cos1(j,k,l) = z1(k)/rad1(j,k,l)
				sin1(j,k,l) = SQRT(x1(j)**2 + y1(k)**2)/rad1(j,k,l)
			END DO
		END DO
   	END DO
END IF

! Get maxmimum distance !
rmax1 = max(maxval(x1), maxval(y1) , maxval(z1))

! Dimensionless distances !
radbar1(:,:,:) = rad1(:,:,:)/rmax1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Get the volume, assuming full box simulation ! 
IF(coordinate_flag == 0) THEN
    vol1(:,:,:) = dh_x * dh_y * dh_z
ELSEIF(coordinate_flag == 1) THEN
   	DO j = -2, nx_1 + 3
    	vol1(j,:,:) = 0.5D0 * (xF1 (j)**2 - xF1 (j-1)**2) * dh_y * dh_z
    ENDDO
ELSEIF(coordinate_flag == 2) THEN
    DO j = -2, nx_1 + 3
	 	DO k = -2, ny_1 + 3
           vol1(j,k,:) = (1.0D0/3.0D0) * (xF1 (j)**3 - xF1 (j-1)**3) * ABS(COS(yF1(k)) - COS(yF1(k-1))) * dh_z
	 	END DO
    ENDDO
END IF

! Special treatment for full box simulation !
IF(NOT(fullx_flag)) THEN
	vol1(:,:,:) = vol1(:,:,:)*2.0D0
END IF
IF(NOT(fully_flag)) THEN
	vol1(:,:,:) = vol1(:,:,:)*2.0D0
END IF
IF(NOT(fullz_flag)) THEN
	vol1(:,:,:) = vol1(:,:,:)*2.0D0
END IF

! Dimensionless volume !
volbar1(:,:,:) = vol1(:,:,:)/rmax1**3

! Set minimum and maximum domain !
nx_min_1 = 1
nx_part_1 = nx_1
ny_min_1 = 1
ny_part_1 = ny_1
nz_min_1 = 1
nz_part_1 = nz_1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE