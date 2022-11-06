!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine prepares the data for spatial discretization,
! due ask the WENO_module to do the reconstruction
! and then combines the results for one Runge-Kutta sub-step
! Prototype developed by Wong Ka Wing in 2010 (or before?)
! Extended by Leung Shing Chi in 2016 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SPATIAL
USE DEFINITION
USE RIEMANN_MODULE
USE WENO_MODULE
USE PPM_MODULE
IMPLICIT NONE

! dummy variables !
INTEGER :: i, j, k, l, p

! Limits of density to be considered
REAL (DP) :: rho_min1, rho_min2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Find frame velocity 
IF (movinggriddm_flag) THEN
	!CALL FINDFRAMEVEL_DM
END IF
IF (movinggridnm_flag) THEN
	!CALL FINDFRAMEVEL_NM
END IF

! Find effective speed !
CALL FINDALPHA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialize source term and fluxes !
IF(RUNDM_flag) THEN
	DO CONCURRENT (j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1, i = imin1:imax1)
		l1(j,k,l,i) = 0.0D0
		sa1(j,k,l,i) = 0.0D0
		sb1(j,k,l,i) = 0.0d0
	END DO
END IF

! For NM !
DO CONCURRENT (j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2, i = imin2:imax2)
	l2(j,k,l,i) = 0.0D0
	sa2(j,k,l,i) = 0.0D0
	sb2(j,k,l,i) = 0.0d0
END DO

! For NM !
DO CONCURRENT (j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2, i = imin2:imax2, p = 1:n_dim)
	flux_2(j,k,l,i,p) = 0.0D0
	dflux_2(j,k,l,i,p) = 0.0d0
END DO

! For DM !
IF(RUNDM_flag) THEN
	DO CONCURRENT (j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1, i = imin1:imax1, p = 1:n_dim)
		flux_1(j,k,l,i,p) = 0.0D0
		dflux_1(j,k,l,i,p) = 0.0d0
	END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section 1: Reconstructions and build states in the x-directions !

! Reconstruct using WENO !
IF(WENO_flag) THEN
	CALL WENO_Reconstruct
END IF

! Build the fluxes and states !
CALL BUILDSTATES

! Choose an appropriate riemann solver !
IF(RUNDM_flag) THEN
	IF(LF_flag) THEN
		CALL LFDM
	END IF
END IF

! For NM !
IF(LF_flag) THEN
	CALL LFNM
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section 2: Get source terms

! For DM sector
IF(runDM_flag) THEN

	! minimum density !
	rho_min1 = 1.1D0 * prim1_a(irho1)

	! Gravitational source term !
	IF(w_gravity) THEN
		DO CONCURRENT (j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1)
			IF(prim1(j,k,l,irho1) > rho_min1) THEN
				sb1(j,k,l,ivel1_x) = prim1(j,k,l,irho1) * dphi_1(j,k,l,x_dir)
				IF(n_dim > 1) THEN
					IF(coordinate_flag == 0)THEN
						sb1(j,k,l,ivel1_y) = prim1(j,k,l,irho1) * dphi_1(j,k,l,y_dir)
					ELSEIF(coordinate_flag == 1)THEN
						sb1(j,k,l,ivel1_y) = prim1(j,k,l,irho1) * dphi_1(j,k,l,y_dir)
					ELSEIF(coordinate_flag == 2)THEN
						sb1(j,k,l,ivel1_y) = prim1(j,k,l,irho1) * dphi_1(j,k,l,y_dir)/x1(j)
					END IF
				ELSEIF(n_dim > 2) THEN
					IF(coordinate_flag == 0)THEN
						sb1(j,k,l,ivel1_z) = prim1(j,k,l,irho1) * dphi_1(j,k,l,z_dir)
					ELSEIF(coordinate_flag == 1)THEN
						sb1(j,k,l,ivel1_z) = prim1(j,k,l,irho1) * dphi_1(j,k,l,z_dir)/x1(j)
					ELSEIF(coordinate_flag == 2)THEN
						sb1(j,k,l,ivel1_z) = prim1(j,k,l,irho1) * dphi_1(j,k,l,z_dir)/(x1(j)*sin1(j,k,l))
					END IF
				END IF
			END IF
		ENDDO
	END IF

	! Add geometric source term to momentum equation according to the coordinate system !
	IF(coordinate_flag == 1) THEN
		DO CONCURRENT (j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1)
			sa1(j,k,l,ivel1_x) = - (p1(j,k,l) + prim1(j,k,l,irho1)*prim1(j,k,l,ivel1_z)**2)
			sa1(j,k,l,ivel1_z) = prim1(j,k,l,irho1)*prim1(j,k,l,ivel1_x)*prim1(j,k,l,ivel1_z)
		END DO
	ELSEIF(coordinate_flag == 2) THEN
		DO CONCURRENT (j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1)
			sa1(j,k,l,ivel1_x) = - (2.0D0*p1(j,k,l) + prim1(j,k,l,irho1)*(prim1(j,k,l,ivel1_y)**2 + prim1(j,k,l,ivel1_z)**2))
			sa1(j,k,l,ivel1_y) = prim1(j,k,l,irho1)*prim1(j,k,l,ivel1_x)*prim1(j,k,l,ivel1_y) 
			sa1(j,k,l,ivel1_z) = prim1(j,k,l,irho1)*prim1(j,k,l,ivel1_x)*prim1(j,k,l,ivel1_z) 
			IF(n_dim > 1) THEN
				sa1(j,k,l,ivel1_y) = sa1(j,k,l,ivel1_y) - (p1(j,k,l) + prim1(j,k,l,irho1)*prim1(j,k,l,ivel1_z)**2)*cos1(j,k,l)/sin1(j,k,l)
				sa1(j,k,l,ivel1_z) = sa1(j,k,l,ivel1_z) + prim1(j,k,l,irho1)*prim1(j,k,l,ivel1_y)*prim1(j,k,l,ivel1_z)*cos1(j,k,l)/sin1(j,k,l)
			END IF
		END DO
	END IF

	! Moving grid !
	IF(movinggriddm_flag) THEN
		DO CONCURRENT (j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1, i = imin1:imax1)
			sb1 (j,k,l,i) = sb1 (j,k,l,i) + cons1 (j,k,l,i) * (3.0D0 * vel1_max / radius1)
		END DO
	END IF

END IF

! Threshold density !
rho_min2 = 1.1D0 * prim2_a(irho2)

! Gravitational source term !
IF(w_gravity) THEN
	DO CONCURRENT (j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2)
		IF(prim2(j,k,l,irho2) > rho_min2) THEN
			sb2(j,k,l,ivel2_x) = prim2(j,k,l,irho2) * dphi_2(j,k,l,x_dir)
			sb2(j,k,l,itau2) = prim2(j,k,l,ivel2_x) * dphi_2(j,k,l,x_dir)
			IF(n_dim > 1) THEN
				IF(coordinate_flag == 0)THEN
					sb2(j,k,l,ivel2_y) = prim2(j,k,l,irho2) * dphi_2(j,k,l,y_dir)
					sb2(j,k,l,itau2) = sb2(j,k,l,itau2) + prim2(j,k,l,ivel2_y) * dphi_2(j,k,l,y_dir)
				ELSEIF(coordinate_flag == 1)THEN
					sb2(j,k,l,ivel2_y) = prim2(j,k,l,irho2) * dphi_2(j,k,l,y_dir)
					sb2(j,k,l,itau2) = sb2(j,k,l,itau2) + prim2(j,k,l,ivel2_y) * dphi_2(j,k,l,y_dir)
				ELSEIF(coordinate_flag == 2)THEN
					sb2(j,k,l,ivel2_y) = prim2(j,k,l,irho2) * dphi_2(j,k,l,y_dir)/x2(j)
					sb2(j,k,l,itau2) = sb2(j,k,l,itau2) + prim2(j,k,l,ivel2_y) * dphi_2(j,k,l,y_dir)/x2(j)
				END IF				
			ELSEIF(n_dim > 2) THEN
				IF(coordinate_flag == 0)THEN
					sb2(j,k,l,ivel2_z) = prim2(j,k,l,irho2) * dphi_2(j,k,l,z_dir)
					sb2(j,k,l,itau2) = sb2(j,k,l,itau2) + prim2(j,k,l,ivel2_z) * dphi_2(j,k,l,z_dir)
				ELSEIF(coordinate_flag == 1)THEN
					sb2(j,k,l,ivel2_z) = prim2(j,k,l,irho2) * dphi_2(j,k,l,z_dir)/x2(j)
					sb2(j,k,l,itau2) = sb2(j,k,l,itau2) + prim2(j,k,l,ivel2_z) * dphi_2(j,k,l,z_dir)/x2(j)
				ELSEIF(coordinate_flag == 2)THEN
					sb2(j,k,l,ivel2_z) = prim2(j,k,l,irho2) * dphi_2(j,k,l,z_dir)/(x2(j)*sin2(j,k,l))
					sb2(j,k,l,itau2) = sb2(j,k,l,itau2) + prim2(j,k,l,ivel2_z) * dphi_2(j,k,l,z_dir)/(x2(j)*sin2(j,k,l))
				END IF
			END IF
			sb2(j,k,l,itau2) = sb2(j,k,l,itau2) * prim2(j,k,l,irho2)
		END IF
	ENDDO
END IF

! Choose coordinate system, add geometric source 
IF(coordinate_flag == 1) THEN
	DO CONCURRENT (j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2)
		sa2(j,k,l,ivel2_x) = - (prim2(j,k,l,itau2) + prim2(j,k,l,irho2)*prim2(j,k,l,ivel2_z)**2)
		sa2(j,k,l,ivel2_z) = prim2(j,k,l,irho2)*prim2(j,k,l,ivel2_x)*prim2(j,k,l,ivel2_z)
	END DO
ELSEIF(coordinate_flag == 2) THEN
	DO CONCURRENT (j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2)
		sa2(j,k,l,ivel2_x) = - (2.0D0*prim2(j,k,l,itau2) + prim2(j,k,l,irho2)*(prim2(j,k,l,ivel2_y)**2 + prim2(j,k,l,ivel2_z)**2))
		sa2(j,k,l,ivel2_y) = prim2(j,k,l,irho2)*prim2(j,k,l,ivel2_x)*prim2(j,k,l,ivel2_y) 
		sa2(j,k,l,ivel2_z) = prim2(j,k,l,irho2)*prim2(j,k,l,ivel2_x)*prim2(j,k,l,ivel2_z) 
		IF(n_dim > 1) THEN
			sa2(j,k,l,ivel2_y) = sa2(j,k,l,ivel2_y) - (prim2(j,k,l,itau2) + prim2(j,k,l,irho2)*prim2(j,k,l,ivel2_z)**2)*cos2(j,k,l)/sin2(j,k,l)
			sa2(j,k,l,ivel2_z) = sa2(j,k,l,ivel2_z) + prim2(j,k,l,irho2)*prim2(j,k,l,ivel2_y)*prim2(j,k,l,ivel2_z)*cos2(j,k,l)/sin2(j,k,l)
		END IF
	END DO
END IF

! dual energy extra source !
IF(dual_energy) THEN

	! Pressure gradient !
	CALL FINDGRADP

	IF(coordinate_flag == 0) THEN
		DO CONCURRENT (j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2)
			sb2(j,k,l,ieps2) = - prim2(j,k,l,ivel2_x)*dp_2 (j,k,l,x_dir)
			IF(n_dim > 1) THEN 
				sb2(j,k,l,ieps2) = sb2(j,k,l,ieps2) - prim2(j,k,l,ivel2_y)*dp_2 (j,k,l,y_dir)
			END IF
			IF(n_dim > 2) THEN 
				sb2(j,k,l,ieps2) = sb2(j,k,l,ieps2) - prim2(j,k,l,ivel2_z)*dp_2 (j,k,l,z_dir)
			END IF
		END DO
	ELSEIF(coordinate_flag == 1) THEN
		DO CONCURRENT (j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2)
			sb2(j,k,l,ieps2) = - prim2(j,k,l,ivel2_x)*dp_2 (j,k,l,x_dir)
			IF(n_dim > 1) THEN 
				sb2(j,k,l,ieps2) = sb2(j,k,l,ieps2) - prim2(j,k,l,ivel2_y)*dp_2 (j,k,l,y_dir)
			END IF
			IF(n_dim > 2) THEN 
				sb2(j,k,l,ieps2) = - prim2(j,k,l,ivel2_z)*dp_2 (j,k,l,z_dir)/x2(j)
			END IF
		END DO
	ELSEIF(coordinate_flag == 1) THEN
		DO CONCURRENT (j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2)
			sb2(j,k,l,ieps2) = - prim2(j,k,l,ivel2_x)*dp_2 (j,k,l,x_dir)
			IF(n_dim > 1) THEN 
			 	sb2(j,k,l,ieps2) = sb2(j,k,l,ieps2) - prim2(j,k,l,ivel2_y)*dp_2 (j,k,l,y_dir)/x2(j)
			END IF
			IF(n_dim > 2) THEN 
			 	sb2(j,k,l,ieps2) = sb2(j,k,l,ieps2) - prim2(j,k,l,ivel2_z)*dp_2 (j,k,l,z_dir)/x2(j)/sin2(j,k,l)
			END IF
		END DO
	END IF
END IF

! Moving grid !
IF(movinggridnm_flag) THEN
	DO CONCURRENT (j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2, i = imin2:imax2)
		sb2 (j,k,l,i) = sb2 (j,k,l,i) + cons2 (j,k,l,i) * (3.0D0 * vel2_max / radius2)
	END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section 4: Get flux derivatives

! Now calculate dfdx accordingly to the corrected flux, do it case by case
IF(coordinate_flag == 0) THEN

	! Single dimension !
	IF(RUNDM_flag) THEN
		DO CONCURRENT (j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1, i = imin1:imax1)
			dflux_1 (j,k,l,i,x_dir) = (flux_1 (j,k,l,i,x_dir) - flux_1 (j-1,k,l,i,x_dir)) / dr1(j)
			IF(n_dim > 1) THEN
				dflux_1 (j,k,l,i,y_dir) = (flux_1 (j,k,l,i,y_dir) - flux_1 (j,k-1,l,i,y_dir)) / dy1
			END IF
			IF(n_dim > 2) THEN
				dflux_1 (j,k,l,i,z_dir) = (flux_1 (j,k,l,i,z_dir) - flux_1 (j,k,l-1,i,z_dir)) / dz1
			END IF
		END DO
	END IF


	! Do for NM !
	DO CONCURRENT (j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2, i = imin2:imax2)
		dflux_2 (j,k,l,i,x_dir) = (flux_2 (j,k,l,i,x_dir) - flux_2 (j-1,k,l,i,x_dir)) / dr2(j)
		IF(n_dim > 1) THEN
			dflux_2 (j,k,l,i,y_dir) = (flux_2 (j,k,l,i,y_dir) - flux_2 (j,k-1,l,i,y_dir)) / dy2
		END IF
		IF(n_dim > 2) THEN
			dflux_2 (j,k,l,i,z_dir) = (flux_2 (j,k,l,i,z_dir) - flux_2 (j,k,l-1,i,z_dir)) / dz2
		END IF
	END DO

ELSEIF (coordinate_flag == 1) THEN

	! along the r-direction !
	IF(RUNDM_flag) THEN
		DO CONCURRENT (j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1, i = imin1:imax1)
			dflux_1 (j,k,l,i,x_dir) = (xF1(j)*flux_1 (j,k,l,i,x_dir) - xF1(j-1)*flux_1 (j-1,k,l,i,x_dir)) / (x1(j)*dr1(j))
			IF(n_dim > 1) THEN
				dflux_1 (j,k,l,i,y_dir) = (flux_1 (j,k,l,i,y_dir) - flux_1 (j,k-1,l,i,y_dir)) / dy1
			END IF
			IF(n_dim > 2) THEN
				dflux_1 (j,k,l,i,z_dir) = (flux_1 (j,k,l,i,z_dir) - flux_1 (j,k,l-1,i,z_dir)) / x1(j)/dz1
			END IF
		END DO
	END IF

	! Do for NM !
	DO CONCURRENT (j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2, i = imin2:imax2)
		dflux_2 (j,k,l,i,x_dir) = (xF2(j)*flux_2 (j,k,l,i,x_dir) - xF2(j-1)*flux_2 (j-1,k,l,i,x_dir)) / (x2(j)*dr2(j))
		IF(n_dim > 1) THEN
			dflux_2 (j,k,l,i,y_dir) = (flux_2 (j,k,l,i,y_dir) - flux_2 (j,k-1,l,i,y_dir)) / dy2
		END IF
		IF(n_dim > 2) THEN
			dflux_2 (j,k,l,i,z_dir) = (flux_2 (j,k,l,i,z_dir) - flux_2 (j,k,l-1,i,z_dir)) / x2(j)/dz2
		END IF
	END DO

ELSEIF (coordinate_flag == 2) THEN

	! along the r-direction !
	IF(RUNDM_flag) THEN
		DO CONCURRENT (j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1, i = imin1:imax1)
			dflux_1 (j,k,l,i,x_dir) = (xF1(j)**2*flux_1 (j,k,l,i,x_dir) - xF1(j-1)**2*flux_1 (j-1,k,l,i,x_dir)) / (x1(j)**2*dr1(j))
			IF(n_dim > 1) THEN
				dflux_1 (j,k,l,i,y_dir) = (DSIN(yF1(k))*flux_1 (j,k,l,i,y_dir) - DSIN(yF1(k-1))*flux_1 (j,k-1,l,i,y_dir)) / (x1(j)*ABS(DCOS(yF1(k)) - DCOS(yF1(k-1))))
			END IF
			IF(n_dim > 2) THEN
				dflux_1 (j,k,l,i,z_dir) = dy1/(ABS(DCOS(yF1(k)) - DCOS(yF1(k-1))))*(flux_1 (j,k,l,i,z_dir) - flux_1 (j,k,l-1,i,z_dir)) / x1(j)/dz1
			END IF
		END DO

	END IF
	
	! Do for NM !
	DO CONCURRENT (j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2, i = imin2:imax2)
		dflux_2 (j,k,l,i,x_dir) = (xF2(j)**2*flux_2 (j,k,l,i,x_dir) - xF2(j-1)**2*flux_2 (j-1,k,l,i,x_dir)) / (x2(j)**2*dr2(j))
		IF(n_dim > 1) THEN
			dflux_2 (j,k,l,i,y_dir) = (DSIN(yF2(k))*flux_2 (j,k,l,i,y_dir) - DSIN(yF2(k-1))*flux_2 (j,k-1,l,i,y_dir)) / (x2(j)*ABS(DCOS(yF2(k)) - DCOS(yF2(k-1))))
		END IF
		IF(n_dim > 2) THEN
			dflux_2 (j,k,l,i,z_dir) = dy2/(ABS(DCOS(yF2(k)) - DCOS(yF2(k-1))))*(flux_2 (j,k,l,i,z_dir) - flux_2 (j,k,l-1,i,z_dir)) / x2(j)/dz2
		END IF
	END DO

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section 5: Calculate l (SSPRK54)

! assign !
IF(coordinate_flag == 0) THEN
	IF(runDM_flag) THEN
		DO CONCURRENT (j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1, i = imin1:imax1)
			l1(j,k,l,i) = - dflux_1 (j,k,l,i,x_dir) - sb1 (j,k,l,i) 
			IF(n_dim > 1) THEN
				l1(j,k,l,i) = l1(j,k,l,i) - dflux_1 (j,k,l,i,y_dir)
			END IF
			IF(n_dim > 2) THEN
				l1(j,k,l,i) = l1(j,k,l,i) - dflux_1 (j,k,l,i,z_dir)
			END IF
		END DO
   	ENDIF
	DO CONCURRENT (j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2, i = imin2:imax2)
   		l2(j,k,l,i) = - dflux_2 (j,k,l,i,x_dir) - sb2 (j,k,l,i)
		IF(n_dim > 1) THEN
			l2(j,k,l,i) = l2(j,k,l,i) - dflux_2 (j,k,l,i,y_dir)
		END IF
		IF(n_dim > 2) THEN
			l2(j,k,l,i) = l2(j,k,l,i) - dflux_2 (j,k,l,i,z_dir)
		END IF
	END DO
ELSE
   	IF(runDM_flag) THEN
        DO CONCURRENT (j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1, i = imin1:imax1)
            l1(j,k,l,i) = - dflux_1 (j,k,l,i,x_dir) - sa1 (j,k,l,i) / x1(j) - sb1 (j,k,l,i)
			IF(n_dim > 1) THEN 
				l1(j,k,l,i) = l1(j,k,l,i) - dflux_1 (j,k,l,i,y_dir)
			END IF
			IF(n_dim > 2) THEN
				l1(j,k,l,i) = l1(j,k,l,i) - dflux_1 (j,k,l,i,z_dir)
			END IF
        ENDDO
   	ENDIF
   	DO CONCURRENT (j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2, i = imin2:imax2)
        l2(j,k,l,i) = - dflux_2 (j,k,l,i,x_dir) - sa2 (j,k,l,i) / x2(j) - sb2 (j,k,l,i)
		IF(n_dim > 1) THEN 
			l2(j,k,l,i) = l2(j,k,l,i) - dflux_2 (j,k,l,i,y_dir)
		END IF
		IF(n_dim > 2) THEN
			l2(j,k,l,i) = l2(j,k,l,i) - dflux_2 (j,k,l,i,z_dir)
		END IF
   	ENDDO
ENDIF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine built the states for left and right edges for the horizontal directions !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BUILDSTATES
USE OMP_LIB
USE DEFINITION 
USE RIEMANN_MODULE
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k, l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Build conserved variables !
If (RUNDM_flag) then
	DO CONCURRENT (j = nx_min_1-1:nx_part_2, k = ny_min_2-1:ny_part_2, l = nz_min_2-1:nz_part_2)
		DO i = imin1, imax1
			uL1 (j,k,l,i,x_dir) = primL1(j,k,l,i,x_dir)*primL1(j,k,l,irho1,x_dir)
			uR1 (j,k,l,i,x_dir) = primR1(j,k,l,i,x_dir)*primR1(j,k,l,irho1,x_dir)
			IF(n_dim > 1) THEN
				uL1 (j,k,l,i,y_dir) = primL1(j,k,l,i,y_dir)*primL1(j,k,l,irho1,y_dir)
				uR1 (j,k,l,i,y_dir) = primR1(j,k,l,i,y_dir)*primR1(j,k,l,irho1,y_dir)
			END IF
			IF(n_dim > 2) THEN
				uL1 (j,k,l,i,z_dir) = primL1(j,k,l,i,z_dir)*primL1(j,k,l,irho1,z_dir)
				uR1 (j,k,l,i,z_dir) = primR1(j,k,l,i,z_dir)*primR1(j,k,l,irho1,z_dir)
			END IF
		END DO
		uL1 (j,k,l,irho1,x_dir) = primL1(j,k,l,irho1,x_dir)
		uR1 (j,k,l,irho1,x_dir) = primR1(j,k,l,irho1,x_dir)
		IF(n_dim > 1) THEN
			uL1 (j,k,l,irho1,y_dir) = primL1(j,k,l,irho1,y_dir)
			uR1 (j,k,l,irho1,y_dir) = primR1(j,k,l,irho1,y_dir)
		END IF
		IF(n_dim > 2) THEN
			uL1 (j,k,l,irho1,z_dir) = primL1(j,k,l,irho1,z_dir)
			uR1 (j,k,l,irho1,z_dir) = primR1(j,k,l,irho1,z_dir)
		END IF

		! Build fluxes For DM !
		DO i = imin1, imax1
			fluxL1 (j,k,l,i,x_dir) = uL1 (j,k,l,i,x_dir) * primL1(j,k,l,ivel1_x,x_dir)
			fluxR1 (j,k,l,i,x_dir) = uR1 (j,k,l,i,x_dir) * primR1(j,k,l,ivel1_x,x_dir)
			IF(n_dim > 1) THEN
				fluxL1 (j,k,l,i,y_dir) = uL1 (j,k,l,i,y_dir) * primL1(j,k,l,ivel1_y,y_dir)
				fluxR1 (j,k,l,i,y_dir) = uR1 (j,k,l,i,y_dir) * primR1(j,k,l,ivel1_y,y_dir)
			END IF
			IF(n_dim > 2) THEN
				fluxL1 (j,k,l,i,z_dir) = uL1 (j,k,l,i,z_dir) * primL1(j,k,l,ivel1_z,z_dir)
				fluxR1 (j,k,l,i,z_dir) = uR1 (j,k,l,i,z_dir) * primR1(j,k,l,ivel1_z,z_dir)
			END IF
		END DO

   		! Add the pressure term to r-momentum equation !
    	fluxL1 (j,k,l,ivel1_x,x_dir) = fluxL1 (j,k,l,ivel1_x,x_dir) + p1L(j,k,l,x_dir)
		fluxR1 (j,k,l,ivel1_x,x_dir) = fluxR1 (j,k,l,ivel1_x,x_dir) + p1R(j,k,l,x_dir)
		IF(n_dim > 1) THEN
    		fluxL1 (j,k,l,ivel1_y,y_dir) = fluxL1 (j,k,l,ivel1_y,y_dir) + p1L(j,k,l,y_dir)
			fluxR1 (j,k,l,ivel1_y,y_dir) = fluxR1 (j,k,l,ivel1_y,y_dir) + p1R(j,k,l,y_dir)
		END IF
		IF(n_dim > 2) THEN
    		fluxL1 (j,k,l,ivel1_z,z_dir) = fluxL1 (j,k,l,ivel1_z,z_dir) + p1L(j,k,l,z_dir)
			fluxR1 (j,k,l,ivel1_z,z_dir) = fluxR1 (j,k,l,ivel1_z,z_dir) + p1R(j,k,l,z_dir)
		END IF

		! Moving grid !
		If(movinggriddm_flag) THEN
         	fluxL1 (j,k,l,i,x_dir) = fluxL1 (j,k,l,i,x_dir) - uL1(j,k,l,i,x_dir)*vf1xL(j)
			fluxR1 (j,k,l,i,x_dir) = fluxR1 (j,k,l,i,x_dir) - uR1(j,k,l,i,x_dir)*vf1xR(j)
			IF(n_dim > 1) THEN
         		fluxL1 (j,k,l,i,y_dir) = fluxL1 (j,k,l,i,y_dir) - uL1(j,k,l,i,y_dir)*vf1yL(k)
				fluxR1 (j,k,l,i,y_dir) = fluxR1 (j,k,l,i,y_dir) - uR1(j,k,l,i,y_dir)*vf1yR(k)
			END IF
			IF(n_dim > 2) THEN
         		fluxL1 (j,k,l,i,z_dir) = fluxL1 (j,k,l,i,z_dir) - uL1(j,k,l,i,z_dir)*vf1zL(l)
				fluxR1 (j,k,l,i,z_dir) = fluxR1 (j,k,l,i,z_dir) - uR1(j,k,l,i,z_dir)*vf1zR(l)
			END IF
		END IF
	END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do the same for NM !
DO CONCURRENT (j = nx_min_2-1:nx_part_2, k = ny_min_2-1:ny_part_2, l = nz_min_2-1:nz_part_2)
	DO i = imin2, imax2
		uL2 (j,k,l,i,x_dir) = primL2(j,k,l,i,x_dir)*primL2(j,k,l,irho2,x_dir)
		uR2 (j,k,l,i,x_dir) = primR2(j,k,l,i,x_dir)*primR2(j,k,l,irho2,x_dir)
		IF(n_dim > 1) THEN
			uL2 (j,k,l,i,y_dir) = primL2(j,k,l,i,y_dir)*primL2(j,k,l,irho2,y_dir)
			uR2 (j,k,l,i,y_dir) = primR2(j,k,l,i,y_dir)*primR2(j,k,l,irho2,y_dir)
		END IF
		IF(n_dim > 2) THEN
			uL2 (j,k,l,i,z_dir) = primL2(j,k,l,i,z_dir)*primL2(j,k,l,irho2,z_dir)
			uR2 (j,k,l,i,z_dir) = primR2(j,k,l,i,z_dir)*primR2(j,k,l,irho2,z_dir)
		END IF
	END DO
	uL2 (j,k,l,irho2,x_dir) = primL2(j,k,l,irho2,x_dir)
	uR2 (j,k,l,irho2,x_dir) = primR2(j,k,l,irho2,x_dir)
	IF(n_dim > 1) THEN
		uL2 (j,k,l,irho2,y_dir) = primL2(j,k,l,irho2,y_dir)
		uR2 (j,k,l,irho2,y_dir) = primR2(j,k,l,irho2,y_dir)
	END IF
	IF(n_dim > 2) THEN
		uL2 (j,k,l,irho2,z_dir) = primL2(j,k,l,irho2,z_dir)
		uR2 (j,k,l,irho2,z_dir) = primR2(j,k,l,irho2,z_dir)
	END IF

	! NM energy equation !
	uL2 (j,k,l,itau2,x_dir) = primL2(j,k,l,irho2,x_dir)*(0.5D0*(primL2(j,k,l,ivel2_x,x_dir)**2 &
						    + primL2(j,k,l,ivel2_y,x_dir)**2 + primL2(j,k,l,ivel2_z,x_dir)**2) + eps2L(j,k,l,x_dir))
	uR2 (j,k,l,itau2,x_dir) = primR2(j,k,l,irho2,x_dir)*(0.5D0*(primR2(j,k,l,ivel2_x,x_dir)**2 &
							+ primR2(j,k,l,ivel2_y,x_dir)**2 + primR2(j,k,l,ivel2_z,x_dir)**2) + eps2R(j,k,l,x_dir))
	IF(n_dim > 1) THEN
		uL2 (j,k,l,itau2,y_dir) = primL2(j,k,l,irho2,y_dir)*(0.5D0*(primL2(j,k,l,ivel2_x,y_dir)**2 &
								+ primL2(j,k,l,ivel2_y,y_dir)**2 + primL2(j,k,l,ivel2_z,y_dir)**2) + eps2L(j,k,l,y_dir))
		uR2 (j,k,l,itau2,y_dir) = primR2(j,k,l,irho2,y_dir)*(0.5D0*(primR2(j,k,l,ivel2_x,y_dir)**2 &
								+ primR2(j,k,l,ivel2_y,y_dir)**2 + primR2(j,k,l,ivel2_z,y_dir)**2) + eps2R(j,k,l,y_dir))
	END IF
	IF(n_dim > 2) THEN
		uL2 (j,k,l,itau2,z_dir) = primL2(j,k,l,irho2,z_dir)*(0.5D0*(primL2(j,k,l,ivel2_x,z_dir)**2 &
								+ primL2(j,k,l,ivel2_y,z_dir)**2 + primL2(j,k,l,ivel2_z,z_dir)**2) + eps2L(j,k,l,z_dir))
		uR2 (j,k,l,itau2,z_dir) = primR2(j,k,l,irho2,z_dir)*(0.5D0*(primR2(j,k,l,ivel2_x,z_dir)**2 &
								+ primR2(j,k,l,ivel2_y,z_dir)**2 + primR2(j,k,l,ivel2_z,z_dir)**2) + eps2R(j,k,l,z_dir))
	END IF

	! dual energy !
	IF(dual_energy) THEN 
		uL2 (j,k,l,ieps2,x_dir) = primL2(j,k,l,ieps2,x_dir)
		uR2 (j,k,l,ieps2,x_dir) = primR2(j,k,l,ieps2,x_dir)
		IF(n_dim > 1) THEN
			uL2 (j,k,l,ieps2,y_dir) = primL2(j,k,l,ieps2,y_dir)
			uR2 (j,k,l,ieps2,y_dir) = primR2(j,k,l,ieps2,y_dir)
		END IF
		IF(n_dim > 2) THEN
			uL2 (j,k,l,ieps2,z_dir) = primL2(j,k,l,ieps2,z_dir)
			uR2 (j,k,l,ieps2,z_dir) = primR2(j,k,l,ieps2,z_dir)
		END IF
	END IF

	! For NM flux !
	DO i = imin2, imax2, 1
		fluxL2 (j,k,l,i,x_dir) = uL2 (j,k,l,i,x_dir) * primL2(j,k,l,ivel2_x,x_dir)
		fluxR2 (j,k,l,i,x_dir) = uR2 (j,k,l,i,x_dir) * primR2(j,k,l,ivel2_x,x_dir)
		IF(n_dim > 1) THEN
			fluxL2 (j,k,l,i,y_dir) = uL2 (j,k,l,i,y_dir) * primL2(j,k,l,ivel2_y,y_dir)
			fluxR2 (j,k,l,i,y_dir) = uR2 (j,k,l,i,y_dir) * primR2(j,k,l,ivel2_y,y_dir)
		END IF
		IF(n_dim > 2) THEN
			fluxL2 (j,k,l,i,z_dir) = uL2 (j,k,l,i,z_dir) * primL2(j,k,l,ivel2_z,z_dir)
			fluxR2 (j,k,l,i,z_dir) = uR2 (j,k,l,i,z_dir) * primR2(j,k,l,ivel2_z,z_dir)
		END IF
	ENDDO

	! Add the pressure term to x-momentum equation !
	fluxL2 (j,k,l,ivel2_x,x_dir) = fluxL2 (j,k,l,ivel2_x,x_dir) + primL2(j,k,l,itau2,x_dir)
	fluxR2 (j,k,l,ivel2_x,x_dir) = fluxR2 (j,k,l,ivel2_x,x_dir) + primR2(j,k,l,itau2,x_dir)
	IF(n_dim > 1) THEN
		fluxL2 (j,k,l,ivel2_y,y_dir) = fluxL2 (j,k,l,ivel2_y,y_dir) + primL2(j,k,l,itau2,y_dir)
		fluxR2 (j,k,l,ivel2_y,y_dir) = fluxR2 (j,k,l,ivel2_y,y_dir) + primR2(j,k,l,itau2,y_dir)
	END IF
	IF(n_dim > 2) THEN
		fluxL2 (j,k,l,ivel2_z,z_dir) = fluxL2 (j,k,l,ivel2_z,z_dir) + primL2(j,k,l,itau2,z_dir)
		fluxR2 (j,k,l,ivel2_z,z_dir) = fluxR2 (j,k,l,ivel2_z,z_dir) + primR2(j,k,l,itau2,z_dir)
	END IF

	! Add the presusre work done term to the energy equation             
	fluxL2 (j,k,l,itau2,x_dir) = fluxL2 (j,k,l,itau2,x_dir) + primL2(j,k,l,itau2,x_dir) * primL2(j,k,l,ivel2_x,x_dir)
	fluxR2 (j,k,l,itau2,x_dir) = fluxR2 (j,k,l,itau2,x_dir) + primR2(j,k,l,itau2,x_dir) * primR2(j,k,l,ivel2_x,x_dir) 
	IF(n_dim > 1) THEN
		fluxL2 (j,k,l,itau2,y_dir) = fluxL2 (j,k,l,itau2,y_dir) + primL2(j,k,l,itau2,y_dir) * primL2(j,k,l,ivel2_y,y_dir)
		fluxR2 (j,k,l,itau2,y_dir) = fluxR2 (j,k,l,itau2,y_dir) + primR2(j,k,l,itau2,y_dir) * primR2(j,k,l,ivel2_y,y_dir) 
	END IF
	IF(n_dim > 2) THEN
		fluxL2 (j,k,l,itau2,z_dir) = fluxL2 (j,k,l,itau2,z_dir) + primL2(j,k,l,itau2,z_dir) * primL2(j,k,l,ivel2_z,z_dir)
		fluxR2 (j,k,l,itau2,z_dir) = fluxR2 (j,k,l,itau2,z_dir) + primR2(j,k,l,itau2,z_dir) * primR2(j,k,l,ivel2_z,z_dir) 
	END IF

	! dual energy !
	IF(dual_energy) THEN  
		fluxL2 (j,k,l,ieps2,x_dir) = fluxL2 (j,k,l,ieps2,x_dir) + primL2(j,k,l,itau2,x_dir) * primL2(j,k,l,ivel2_x,x_dir)
		fluxR2 (j,k,l,ieps2,x_dir) = fluxR2 (j,k,l,ieps2,x_dir) + primR2(j,k,l,itau2,x_dir) * primR2(j,k,l,ivel2_x,x_dir)
		IF(n_dim > 1) THEN
			fluxL2 (j,k,l,ieps2,y_dir) = fluxL2 (j,k,l,ieps2,y_dir) + primL2(j,k,l,itau2,y_dir) * primL2(j,k,l,ivel2_y,y_dir)
			fluxR2 (j,k,l,ieps2,y_dir) = fluxR2 (j,k,l,ieps2,y_dir) + primR2(j,k,l,itau2,y_dir) * primR2(j,k,l,ivel2_y,y_dir)
		END IF
		IF(n_dim > 2) THEN
			fluxL2 (j,k,l,ieps2,z_dir) = fluxL2 (j,k,l,ieps2,z_dir) + primL2(j,k,l,itau2,z_dir) * primL2(j,k,l,ivel2_z,z_dir)
			fluxR2 (j,k,l,ieps2,z_dir) = fluxR2 (j,k,l,ieps2,z_dir) + primR2(j,k,l,itau2,z_dir) * primR2(j,k,l,ivel2_z,z_dir)
		END IF
	END IF

	! Extra flux term for moving grid !
	If(movinggridnm_flag) THEN
        fluxL2 (j,k,l,i,x_dir) = fluxL2 (j,k,l,i,x_dir) - uL2(j,k,l,i,x_dir)*vf2xL(j)
		fluxR2 (j,k,l,i,x_dir) = fluxR2 (j,k,l,i,x_dir) - uR2(j,k,l,i,x_dir)*vf2xR(j)
		IF(n_dim > 1) THEN
       	 	fluxL2 (j,k,l,i,y_dir) = fluxL2 (j,k,l,i,y_dir) - uL2(j,k,l,i,y_dir)*vf2yL(k)
			fluxR2 (j,k,l,i,y_dir) = fluxR2 (j,k,l,i,y_dir) - uR2(j,k,l,i,y_dir)*vf2yR(k)
		END IF
		IF(n_dim > 2) THEN
       	 	fluxL2 (j,k,l,i,z_dir) = fluxL2 (j,k,l,i,z_dir) - uL2(j,k,l,i,z_dir)*vf2zL(l)
			fluxR2 (j,k,l,i,z_dir) = fluxR2 (j,k,l,i,z_dir) - uR2(j,k,l,i,z_dir)*vf2zR(l)
		END IF
	END IF
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE