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

SUBROUTINE GetGridNM
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
        x2(j) = (DBLE(k - n_mid) - 0.5D0) * dh_x
		xF2(j) = DBLE(k - n_mid)* dh_x
      	dr2(j) = dh_x
   	ENDDO
END IF
IF(fully_flag) THEN
	n_mid = ny_2/2
   	DO j = -2, ny_2 + 3
        y2(j) = (DBLE(k - n_mid) - 0.5D0) * dh_y
		yF2(j) = DBLE(k - n_mid)* dh_y
   	ENDDO
END IF
IF(fullz_flag) THEN
	n_mid = nz_2/2
   	DO j = -2, nz_2 + 3
        z2(j) = (DBLE(k - n_mid) - 0.5D0) * dh_z
		zF2(j) = DBLE(k - n_mid)* dh_z
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
 			cos2(j,k,:) = z2(k)/rad2(j,k,:)
			sin2(j,k,:) = r2(k)/rad2(j,k,:)
		END DO
   	END DO
ELSE
   	DO j = -2, nx_2 + 3
		DO k = -2, ny_2 + 3
			DO l = -2, nz_2 + 3
				rad2(j,k,l) = SQRT(x2(j)**2 + y2(k)**2 + z2(k)**2)
 				cos2(j,k,l) = z2(k)/rad2(j,k,l)
				sin2(j,k,l) = r2(k)/rad2(j,k,l)
			END DO
		END DO
   	END DO
END IF

! Get maxmimum distance !
rmax2 = max(maxval(x2), maxval(y2) , maxval(z2))

! Dimensionless distances !
radbar2(:,:,:) = rad2(:,:,:)/rmax2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Get the volume 
IF(coordinate_flag == 0) THEN
    vol2(:,:,:) = 8.0D0 * dh_x * dh_y * dh_z
ELSEIF(coordinate_flag == 1) THEN
   	IF(hemisphere_flag == 1) THEN
      	! You need to change this volume if we use
      	! stricktly Cartesian coordinate
      	DO j = -4, length_step_r_2 + 5, 1
         	vol2(j,:) = pi * (rF2 (j)**2 - rF2 (j-1)**2) * dxz
      	ENDDO
   	ELSEIF(hemisphere_flag == 0) THEN
      	DO j = -4, length_step_r_2 + 5, 1
         	vol2(j,:) = 2.0D0 * pi * (rF2 (j)**2 - rF2 (j-1)**2) * dxz
      	ENDDO
   	ELSE
      	STOP 'Check hemisphere flag, stopped at GetGrid'
   	ENDIF
ELSEIF(coordinate_flag == 2) THEN
    DO j = -4, length_step_r_2 + 5, 1	
	 	DO k = -4, length_step_z_2 + 5, 1
           vol2(j,k) = (4.0D0/3.0D0) * pi * (rF2 (j)**3 - rF2 (j-1)**3) * ABS(COS(zF2(k)) - COS(zF2(k-1)))
	 	END DO
    ENDDO
END IF

! Dimensionless volume !
volbar2(:,:) = vol2(:,:)/rmax2**3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do everything the same, but for DM !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE GetGridDM
USE definition
IMPLICIT NONE

! Dummy variables
INTEGER :: j, k, x

! Real variables !
REAL (DP) :: dxr, dxz

! Assign !
IF(movinggriddm_flag == 1) THEN
	dxr = dx1
	dxz = dx1
ELSE
	dxr = dx1_r
	dxz = dx1_z
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For DM motherboard !

! Get the z-position
IF(coordinate_flag == 0) THEN
	DO k = -4, length_step_z_1 + 5, 1
      	z1(k) = (DBLE(k) - 0.5D0) * dxz
		zF1(k) = DBLE(k)* dxz
   	ENDDO
ELSEIF(coordinate_flag == 1) THEN
   	IF(hemisphere_flag == 1) THEN
      	x = length_step_z_1/2
      	DO k = -4, length_step_z_1 + 5, 1
         	z1(k) = (DBLE(k - x) - 0.5D0) * dxz
			zF1(k) = DBLE(k - x)* dxz
      	ENDDO
   	ELSE
		DO k = -4, length_step_z_1 + 5, 1
         	z1(k) = (DBLE(k) - 0.5D0) * dxz
			zF1(k) = DBLE(k)* dxz
      	ENDDO
   	ENDIF
ELSEIF(coordinate_flag == 2) THEN
	DO k = -4, length_step_z_1 + 5, 1
      	z1(k) = (DBLE(k) - 0.5D0) * dxz
		zF1(k) = DBLE(k)* dxz
   	ENDDO
ENDIF

! Get the r-position
IF(coordinate_flag == 0) THEN
   	DO j = -4, length_step_r_1 + 5, 1
      	r1(j) = (DBLE(j) - 0.5D0) * dxr
		rF1(j) = DBLE(j) * dxr
		dr1(j) = dxr
   	ENDDO
ELSEIF(coordinate_flag == 1) THEN
  	DO j = -4, length_step_r_1 + 5, 1
      	r1(j) = (DBLE(j) - 0.5D0) * dxr
        rF1(j) = DBLE(j) * dxr
		dr1(j) = dxr
   	ENDDO
ELSEIF(coordinate_flag == 2) THEN
   	IF(fornax_flag == 1) THEN
      	DO j = -4, length_step_r_1 + 5, 1
			rF1(j) = A_fornax * xt * dsinh(dble(j)/xt)
      	END DO
      	DO j = -4, length_step_r_1 + 5
        	r1(j) = 0.5D0*(rF1(j) + rF1(j-1))
        	dr1(j) = rF1(j) - rF1(j-1)
      	END DO
   	ELSE
      	DO j = -4, length_step_r_1 + 5, 1
         	r1(j) = (DBLE(j) - 0.5D0) * dxr
         	rF1(j) = DBLE(j) * dxr
         	dr1(j) = dxr
      	ENDDO
   	END IF
ENDIF

! Radial distances !
IF(coordinate_flag == 2) THEN
   	DO j = -4, length_step_r_1 + 5, 1
		DO k = -4, length_step_z_1 + 5, 1
			rad1(j,k) = r1(j)
 			cos1(j,k) = DCOS(z1(k))
			sin1(j,k) = DSIN(z1(k))
		END DO
   	END DO
ELSE
  	DO j = -4, length_step_r_1 + 5, 1
		DO k = -4, length_step_z_1 + 5, 1
			rad1(j,k) = SQRT(r1(j)**2 + z1(k)**2)
			cos1(j,k) = z1(k)/rad1(j,k)
			sin1(j,k) = r1(k)/rad1(j,k)
		END DO
  	END DO
END IF

! Get maxmimum distance !
rmax1 = max(maxval(r1), maxval(z1))

! Dimensionless distances !
radbar1(:,:) = rad1(:,:)/rmax1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For DM !

IF(coordinate_flag == 0) THEN
   	DO j = -4, length_step_r_1 + 5, 1
      	vol1(j,:) = dxr * dxz
   	ENDDO
ELSEIF(coordinate_flag == 1) THEN
   	IF(hemisphere_flag == 1) THEN
      	! You need to change this volume if we use
      	! stricktly Cartesian coordinate
      	DO j = -4, length_step_r_1 + 5, 1
         	vol1(j,:) = pi * (rF1 (j)**2 - rF1 (j-1)**2) * dxz
      	ENDDO
   	ELSEIF(hemisphere_flag == 0) THEN
      	DO j = -4, length_step_r_1 + 5, 1
         	vol1(j,:) = 2.0D0 * pi * (rF1 (j)**2 - rF1 (j-1)**2) * dxz
      	ENDDO
   	ELSE
      	STOP 'Check hemisphere flag, stopped at GetGrid'
   	ENDIF
ELSEIF(coordinate_flag == 2) THEN
    DO j = -4, length_step_r_1 + 5, 1	
	 	DO k = -4, length_step_z_1 + 5, 1
        	vol1(j,k) = (4.0D0/3.0D0) * pi * (rF1 (j)**3 - rF1 (j-1)**3) * ABS(COS(zF1(k)) - COS(zF1(k-1)))
	 	END DO
    ENDDO
END IF

! Dimensionless volume !
volbar1(:,:) = vol1(:,:)/rmax1**3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE