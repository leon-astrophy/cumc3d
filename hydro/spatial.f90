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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! dummy variables !
INTEGER :: i, j, k, l, p

! Limits of density to be considered
REAL (DP) :: rho_min1, rho_min2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section 1: Reconstructions and build states !

! Find effective speed !
CALL FINDALPHA

! Reconstruct using WENO !
IF(WENO_flag) THEN
	CALL WENO_Reconstruct
END IF

! Build the fluxes and states !
CALL BUILDSTATES

! Choose an appropriate riemann solver !
IF(DM_flag) THEN
	IF(LF_flag) THEN
		CALL LFDM
	END IF
END IF

! For NM !
IF(LF_flag) THEN
	CALL LFNM
END IF

! Pressure gradient !
IF(dual_energy) THEN
	CALL FINDGRADP
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section 2: Get source terms, fluxes, and flux graidents for NM

! Threshold density !
rho_min2 = 1.1D0 * prim2_a(irho2)

!  We loop over the entire domain !
DO CONCURRENT (j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2)

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Check for gravitational source terms !
	IF(w_gravity) THEN
		IF(prim2(irho2,j,k,l) > rho_min2) THEN
			sb2(ivel2_x,j,k,l) = prim2(irho2,j,k,l) * dphi_2(x_dir,j,k,l)
			sb2(itau2,j,k,l) = prim2(ivel2_x,j,k,l) * dphi_2(x_dir,j,k,l)
			IF(n_dim > 1) THEN
				IF(coordinate_flag == 0)THEN
					sb2(ivel2_y,j,k,l) = prim2(irho2,j,k,l) * dphi_2(y_dir,j,k,l)
					sb2(itau2,j,k,l) = sb2(itau2,j,k,l) + prim2(ivel2_y,j,k,l) * dphi_2(y_dir,j,k,l)
				ELSEIF(coordinate_flag == 1)THEN
					sb2(ivel2_y,j,k,l) = prim2(irho2,j,k,l) * dphi_2(y_dir,j,k,l)
					sb2(itau2,j,k,l) = sb2(itau2,j,k,l) + prim2(ivel2_y,j,k,l) * dphi_2(y_dir,j,k,l)
				ELSEIF(coordinate_flag == 2)THEN
					sb2(ivel2_y,j,k,l) = prim2(irho2,j,k,l) * dphi_2(y_dir,j,k,l)/x2(j)
					sb2(itau2,j,k,l) = sb2(itau2,j,k,l) + prim2(ivel2_y,j,k,l) * dphi_2(y_dir,j,k,l)/x2(j)
				END IF				
			ELSEIF(n_dim > 2) THEN
				IF(coordinate_flag == 0)THEN
					sb2(ivel2_z,j,k,l) = prim2(irho2,j,k,l) * dphi_2(z_dir,j,k,l)
					sb2(itau2,j,k,l) = sb2(itau2,j,k,l) + prim2(ivel2_z,j,k,l) * dphi_2(z_dir,j,k,l)
				ELSEIF(coordinate_flag == 1)THEN
					sb2(ivel2_z,j,k,l) = prim2(irho2,j,k,l) * dphi_2(z_dir,j,k,l)/x2(j)
					sb2(itau2,j,k,l) = sb2(itau2,j,k,l) + prim2(ivel2_z,j,k,l) * dphi_2(z_dir,j,k,l)/x2(j)
				ELSEIF(coordinate_flag == 2)THEN
					sb2(ivel2_z,j,k,l) = prim2(irho2,j,k,l) * dphi_2(z_dir,j,k,l)/(x2(j)*sin2(j,k,l))
					sb2(itau2,j,k,l) = sb2(itau2,j,k,l) + prim2(ivel2_z,j,k,l) * dphi_2(z_dir,j,k,l)/(x2(j)*sin2(j,k,l))
				END IF
			END IF
			sb2(itau2,j,k,l) = sb2(itau2,j,k,l) * prim2(irho2,j,k,l)
		END IF
	END IF

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Choose coordinate system, add geometric source term !
	IF(coordinate_flag == 1) THEN
		sa2(ivel2_x,j,k,l) = - (prim2(itau2,j,k,l) + prim2(irho2,j,k,l)*prim2(ivel2_z,j,k,l)**2)
		sa2(ivel2_z,j,k,l) = prim2(irho2,j,k,l)*prim2(ivel2_x,j,k,l)*prim2(ivel2_z,j,k,l)
	ELSEIF(coordinate_flag == 2) THEN
		sa2(ivel2_x,j,k,l) = - (2.0D0*prim2(itau2,j,k,l) + prim2(irho2,j,k,l)*(prim2(ivel2_y,j,k,l)**2 + prim2(ivel2_z,j,k,l)**2))
		sa2(ivel2_y,j,k,l) = prim2(irho2,j,k,l)*prim2(ivel2_x,j,k,l)*prim2(ivel2_y,j,k,l) 
		sa2(ivel2_z,j,k,l) = prim2(irho2,j,k,l)*prim2(ivel2_x,j,k,l)*prim2(ivel2_z,j,k,l) 
		IF(n_dim > 1) THEN
			sa2(ivel2_y,j,k,l) = sa2(ivel2_y,j,k,l) - (prim2(itau2,j,k,l) + prim2(irho2,j,k,l)*prim2(ivel2_z,j,k,l)**2)*cos2(j,k,l)/sin2(j,k,l)
			sa2(ivel2_z,j,k,l) = sa2(ivel2_z,j,k,l) + prim2(irho2,j,k,l)*prim2(ivel2_y,j,k,l)*prim2(ivel2_z,j,k,l)*cos2(j,k,l)/sin2(j,k,l)
		END IF
	END IF

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! dual energy extra source !
	IF(dual_energy) THEN

		! Choose according to coordinate system !
		IF(coordinate_flag == 0) THEN
			sb2(ieps2,j,k,l) = - prim2(ivel2_x,j,k,l)*dp_2 (x_dir,j,k,l)
			IF(n_dim > 1) THEN 
				sb2(ieps2,j,k,l) = sb2(ieps2,j,k,l) - prim2(ivel2_y,j,k,l)*dp_2 (y_dir,j,k,l)
			END IF
			IF(n_dim > 2) THEN 
				sb2(ieps2,j,k,l) = sb2(ieps2,j,k,l) - prim2(ivel2_z,j,k,l)*dp_2 (z_dir,j,k,l)
			END IF
		ELSEIF(coordinate_flag == 1) THEN
			sb2(ieps2,j,k,l) = - prim2(ivel2_x,j,k,l)*dp_2 (x_dir,j,k,l)
			IF(n_dim > 1) THEN 
				sb2(ieps2,j,k,l) = sb2(ieps2,j,k,l) - prim2(ivel2_y,j,k,l)*dp_2 (y_dir,j,k,l)
			END IF
			IF(n_dim > 2) THEN 
				sb2(ieps2,j,k,l) = - prim2(ivel2_z,j,k,l)*dp_2 (z_dir,j,k,l)/x2(j)
			END IF
		ELSEIF(coordinate_flag == 2) THEN
			sb2(ieps2,j,k,l) = - prim2(ivel2_x,j,k,l)*dp_2 (x_dir,j,k,l)
			IF(n_dim > 1) THEN 
			 	sb2(ieps2,j,k,l) = sb2(ieps2,j,k,l) - prim2(ivel2_y,j,k,l)*dp_2 (y_dir,j,k,l)/x2(j)
			END IF
			IF(n_dim > 2) THEN 
			 	sb2(ieps2,j,k,l) = sb2(ieps2,j,k,l) - prim2(ivel2_z,j,k,l)*dp_2 (z_dir,j,k,l)/x2(j)/sin2(j,k,l)
			END IF
		END IF

	END IF

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Now calculate flux gradient, choose according to the coordinate system !
	IF(coordinate_flag == 0) THEN
		DO i = imin2, imax2
			dflux_2 (i,x_dir,j,k,l) = (flux_2 (i,x_dir,j,k,l) - flux_2 (i,x_dir,j-1,k,l)) / dr2(j)
			IF(n_dim > 1) THEN
				dflux_2 (i,y_dir,j,k,l) = (flux_2 (i,y_dir,j,k,l) - flux_2 (i,y_dir,j,k-1,l)) / dy2
			END IF
			IF(n_dim > 2) THEN
				dflux_2 (i,z_dir,j,k,l) = (flux_2 (i,z_dir,j,k,l) - flux_2 (i,z_dir,j,k,l-1)) / dz2
			END IF
		END DO
	ELSEIF(coordinate_flag == 1) THEN
		DO i = imin2, imax2
			dflux_2 (i,x_dir,j,k,l) = (xF2(j)*flux_2 (i,x_dir,j,k,l) - xF2(j-1)*flux_2 (i,x_dir,j-1,k,l)) / (x2(j)*dr2(j))
			IF(n_dim > 1) THEN
				dflux_2 (i,y_dir,j,k,l) = (flux_2 (i,y_dir,j,k,l) - flux_2 (i,y_dir,j,k-1,l)) / dy2
			END IF
			IF(n_dim > 2) THEN
				dflux_2 (i,z_dir,j,k,l) = (flux_2 (i,z_dir,j,k,l) - flux_2 (i,z_dir,j,k,l-1)) / x2(j)/dz2
			END IF
		END DO
	ELSEIF(coordinate_flag == 2) THEN
		DO i = imin2, imax2
			dflux_2 (i,x_dir,j,k,l) = (xF2(j)**2*flux_2 (i,x_dir,j,k,l) - xF2(j-1)**2*flux_2 (i,x_dir,j-1,k,l)) / (x2(j)**2*dr2(j))
			IF(n_dim > 1) THEN
				dflux_2 (i,y_dir,j,k,l) = (DSIN(yF2(k))*flux_2 (i,y_dir,j,k,l) - DSIN(yF2(k-1))*flux_2 (i,y_dir,j,k-1,l)) / (x2(j)*ABS(DCOS(yF2(k)) - DCOS(yF2(k-1))))
			END IF
			IF(n_dim > 2) THEN
				dflux_2 (i,z_dir,j,k,l) = dy2/(ABS(DCOS(yF2(k)) - DCOS(yF2(k-1))))*(flux_2 (i,z_dir,j,k,l) - flux_2 (i,z_dir,j,k,l-1)) / x2(j)/dz2
			END IF
		END DO
	END IF

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Final step, get rungekutta operator, LHS of the hydro equation !
	IF(coordinate_flag == 0) THEN
		DO i = imin2, imax2
   			l2(i,j,k,l) = - dflux_2 (i,x_dir,j,k,l) - sb2 (i,j,k,l)
			IF(n_dim > 1) THEN
				l2(i,j,k,l) = l2(i,j,k,l) - dflux_2 (i,y_dir,j,k,l)
			END IF
			IF(n_dim > 2) THEN
				l2(i,j,k,l) = l2(i,j,k,l) - dflux_2 (i,z_dir,j,k,l)
			END IF
		END DO
   	ELSE
		DO i = imin2, imax2
        	l2(i,j,k,l) = - dflux_2 (i,x_dir,j,k,l) - sa2 (i,j,k,l) / x2(j) - sb2 (i,j,k,l)
			IF(n_dim > 1) THEN 
				l2(i,j,k,l) = l2(i,j,k,l) - dflux_2 (i,y_dir,j,k,l)
			END IF
			IF(n_dim > 2) THEN
				l2(i,j,k,l) = l2(i,j,k,l) - dflux_2 (i,z_dir,j,k,l)
			END IF
		END DO
   	END IF

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section 3: repeat everything for DM !

! For DM sector
IF(DM_flag) THEN

	! minimum density !
	rho_min1 = 1.1D0 * prim1_a(irho1)

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Gravitational source term !
	DO CONCURRENT (j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1)
		IF(w_gravity) THEN
			IF(prim1(irho1,j,k,l) > rho_min1) THEN
				sb1(ivel1_x,j,k,l) = prim1(irho1,j,k,l) * dphi_1(x_dir,j,k,l)
				IF(n_dim > 1) THEN
					IF(coordinate_flag == 0)THEN
						sb1(ivel1_y,j,k,l) = prim1(irho1,j,k,l) * dphi_1(y_dir,j,k,l)
					ELSEIF(coordinate_flag == 1)THEN
						sb1(ivel1_y,j,k,l) = prim1(irho1,j,k,l) * dphi_1(y_dir,j,k,l)
					ELSEIF(coordinate_flag == 2)THEN
						sb1(ivel1_y,j,k,l) = prim1(irho1,j,k,l) * dphi_1(y_dir,j,k,l)/x1(j)
					END IF
				ELSEIF(n_dim > 2) THEN
					IF(coordinate_flag == 0)THEN
						sb1(ivel1_z,j,k,l) = prim1(irho1,j,k,l) * dphi_1(z_dir,j,k,l)
					ELSEIF(coordinate_flag == 1)THEN
						sb1(ivel1_z,j,k,l) = prim1(irho1,j,k,l) * dphi_1(z_dir,j,k,l)/x1(j)
					ELSEIF(coordinate_flag == 2)THEN
						sb1(ivel1_z,j,k,l) = prim1(irho1,j,k,l) * dphi_1(z_dir,j,k,l)/(x1(j)*sin1(j,k,l))
					END IF
				END IF
			END IF
		END IF

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Add geometric source term to momentum equation according to the coordinate system !
		IF(coordinate_flag == 1) THEN
			sa1(ivel1_x,j,k,l) = - (p1(j,k,l) + prim1(irho1,j,k,l)*prim1(ivel1_z,j,k,l)**2)
			sa1(ivel1_z,j,k,l) = prim1(irho1,j,k,l)*prim1(ivel1_x,j,k,l)*prim1(ivel1_z,j,k,l)
		ELSEIF(coordinate_flag == 2) THEN
			sa1(ivel1_x,j,k,l) = - (2.0D0*p1(j,k,l) + prim1(irho1,j,k,l)*(prim1(ivel1_y,j,k,l)**2 + prim1(ivel1_z,j,k,l)**2))
			sa1(ivel1_y,j,k,l) = prim1(irho1,j,k,l)*prim1(ivel1_x,j,k,l)*prim1(ivel1_y,j,k,l) 
			sa1(ivel1_z,j,k,l) = prim1(irho1,j,k,l)*prim1(ivel1_x,j,k,l)*prim1(ivel1_z,j,k,l) 
			IF(n_dim > 1) THEN
				sa1(ivel1_y,j,k,l) = sa1(ivel1_y,j,k,l) - (p1(j,k,l) + prim1(irho1,j,k,l)*prim1(ivel1_z,j,k,l)**2)*cos1(j,k,l)/sin1(j,k,l)
				sa1(ivel1_z,j,k,l) = sa1(ivel1_z,j,k,l) + prim1(irho1,j,k,l)*prim1(ivel1_y,j,k,l)*prim1(ivel1_z,j,k,l)*cos1(j,k,l)/sin1(j,k,l)
			END IF
		END IF

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Now calculate dfdx accordingly to the corrected flux, do it case by case
		IF(coordinate_flag == 0) THEN
			DO i = imin1, imax1
				dflux_1 (i,x_dir,j,k,l) = (flux_1 (i,x_dir,j,k,l) - flux_1 (i,x_dir,j-1,k,l)) / dr1(j)
				IF(n_dim > 1) THEN
					dflux_1 (i,y_dir,j,k,l) = (flux_1 (i,y_dir,j,k,l) - flux_1 (i,y_dir,j,k-1,l)) / dy1
				END IF
				IF(n_dim > 2) THEN
					dflux_1 (i,z_dir,j,k,l) = (flux_1 (i,z_dir,j,k,l) - flux_1 (i,z_dir,j,k,l-1)) / dz1
				END IF
			END DO
		ELSEIF (coordinate_flag == 1) THEN
			DO i = imin1, imax1
				dflux_1 (i,x_dir,j,k,l) = (xF1(j)*flux_1 (i,x_dir,j,k,l) - xF1(j-1)*flux_1 (i,x_dir,j-1,k,l)) / (x1(j)*dr1(j))
				IF(n_dim > 1) THEN
					dflux_1 (i,y_dir,j,k,l) = (flux_1 (i,y_dir,j,k,l) - flux_1 (i,y_dir,j,k-1,l)) / dy1
				END IF
				IF(n_dim > 2) THEN
					dflux_1 (i,z_dir,j,k,l) = (flux_1 (i,z_dir,j,k,l) - flux_1 (i,z_dir,j,k,l-1)) / x1(j)/dz1
				END IF
			END DO
		ELSEIF (coordinate_flag == 2) THEN
			DO i = imin1, imax1
				dflux_1 (i,x_dir,j,k,l) = (xF1(j)**2*flux_1 (i,x_dir,j,k,l) - xF1(j-1)**2*flux_1 (i,x_dir,j-1,k,l)) / (x1(j)**2*dr1(j))
				IF(n_dim > 1) THEN
					dflux_1 (i,y_dir,j,k,l) = (DSIN(yF1(k))*flux_1 (i,y_dir,j,k,l) - DSIN(yF1(k-1))*flux_1 (i,y_dir,j,k-1,l)) / (x1(j)*ABS(DCOS(yF1(k)) - DCOS(yF1(k-1))))
				END IF
				IF(n_dim > 2) THEN
					dflux_1 (i,z_dir,j,k,l) = dy1/(ABS(DCOS(yF1(k)) - DCOS(yF1(k-1))))*(flux_1 (i,z_dir,j,k,l) - flux_1 (i,z_dir,j,k,l-1)) / x1(j)/dz1
				END IF
			END DO
		END IF

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Last step, get LHS of the rungekutta operator
		IF(coordinate_flag == 0) THEN
			DO i = imin1, imax1
				l1(i,j,k,l) = - dflux_1 (i,x_dir,j,k,l) - sb1 (i,j,k,l) 
				IF(n_dim > 1) THEN
					l1(i,j,k,l) = l1(i,j,k,l) - dflux_1 (i,y_dir,j,k,l)
				END IF
				IF(n_dim > 2) THEN
					l1(i,j,k,l) = l1(i,j,k,l) - dflux_1 (i,z_dir,j,k,l)
				END IF
			END DO
		ELSE
			DO i = imin1, imax1
            	l1(i,j,k,l) = - dflux_1 (i,x_dir,j,k,l) - sa1 (i,j,k,l) / x1(j) - sb1 (i,j,k,l)
				IF(n_dim > 1) THEN 
					l1(i,j,k,l) = l1(i,j,k,l) - dflux_1 (i,y_dir,j,k,l)
				END IF
				IF(n_dim > 2) THEN
					l1(i,j,k,l) = l1(i,j,k,l) - dflux_1 (i,z_dir,j,k,l)
				END IF
			END DO       
   		ENDIF
	END DO
ENDIF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine built the states for left and right edges for the horizontal directions !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BUILDSTATES
USE DEFINITION 
USE RIEMANN_MODULE
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k, l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Build conserved variables !
If (DM_flag) then
	DO CONCURRENT (j = nx_min_1-1:nx_part_2, k = ny_min_2-1:ny_part_2, l = nz_min_2-1:nz_part_2)
		DO i = imin1, imax1
			uL1 (i,x_dir,j,k,l) = primL1(i,x_dir,j,k,l)*primL1(irho1,x_dir,j,k,l)
			uR1 (i,x_dir,j,k,l) = primR1(i,x_dir,j,k,l)*primR1(irho1,x_dir,j,k,l)
			IF(n_dim > 1) THEN
				uL1 (i,y_dir,j,k,l) = primL1(i,y_dir,j,k,l)*primL1(irho1,y_dir,j,k,l)
				uR1 (i,y_dir,j,k,l) = primR1(i,y_dir,j,k,l)*primR1(irho1,y_dir,j,k,l)
			END IF
			IF(n_dim > 2) THEN
				uL1 (i,z_dir,j,k,l) = primL1(i,z_dir,j,k,l)*primL1(irho1,z_dir,j,k,l)
				uR1 (i,z_dir,j,k,l) = primR1(i,z_dir,j,k,l)*primR1(irho1,z_dir,j,k,l)
			END IF
		END DO
		uL1 (irho1,x_dir,j,k,l) = primL1(irho1,x_dir,j,k,l)
		uR1 (irho1,x_dir,j,k,l) = primR1(irho1,x_dir,j,k,l)
		IF(n_dim > 1) THEN
			uL1 (irho1,y_dir,j,k,l) = primL1(irho1,y_dir,j,k,l)
			uR1 (irho1,y_dir,j,k,l) = primR1(irho1,y_dir,j,k,l)
		END IF
		IF(n_dim > 2) THEN
			uL1 (irho1,z_dir,j,k,l) = primL1(irho1,z_dir,j,k,l)
			uR1 (irho1,z_dir,j,k,l) = primR1(irho1,z_dir,j,k,l)
		END IF

		! Build fluxes For DM !
		DO i = imin1, imax1
			fluxL1 (i,x_dir,j,k,l) = uL1 (i,x_dir,j,k,l) * primL1(ivel1_x,x_dir,j,k,l)
			fluxR1 (i,x_dir,j,k,l) = uR1 (i,x_dir,j,k,l) * primR1(ivel1_x,x_dir,j,k,l)
			IF(n_dim > 1) THEN
				fluxL1 (i,y_dir,j,k,l) = uL1 (i,y_dir,j,k,l) * primL1(ivel1_y,y_dir,j,k,l)
				fluxR1 (i,y_dir,j,k,l) = uR1 (i,y_dir,j,k,l) * primR1(ivel1_y,y_dir,j,k,l)
			END IF
			IF(n_dim > 2) THEN
				fluxL1 (i,z_dir,j,k,l) = uL1 (i,z_dir,j,k,l) * primL1(ivel1_z,z_dir,j,k,l)
				fluxR1 (i,z_dir,j,k,l) = uR1 (i,z_dir,j,k,l) * primR1(ivel1_z,z_dir,j,k,l)
			END IF
		END DO

   		! Add the pressure term to r-momentum equation !
    	fluxL1 (ivel1_x,x_dir,j,k,l) = fluxL1 (ivel1_x,x_dir,j,k,l) + p1L(x_dir,j,k,l)
		fluxR1 (ivel1_x,x_dir,j,k,l) = fluxR1 (ivel1_x,x_dir,j,k,l) + p1R(x_dir,j,k,l)
		IF(n_dim > 1) THEN
    		fluxL1 (ivel1_y,y_dir,j,k,l) = fluxL1 (ivel1_y,y_dir,j,k,l) + p1L(y_dir,j,k,l)
			fluxR1 (ivel1_y,y_dir,j,k,l) = fluxR1 (ivel1_y,y_dir,j,k,l) + p1R(y_dir,j,k,l)
		END IF
		IF(n_dim > 2) THEN
    		fluxL1 (ivel1_z,z_dir,j,k,l) = fluxL1 (ivel1_z,z_dir,j,k,l) + p1L(z_dir,j,k,l)
			fluxR1 (ivel1_z,z_dir,j,k,l) = fluxR1 (ivel1_z,z_dir,j,k,l) + p1R(z_dir,j,k,l)
		END IF
	END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do the same for NM !
DO CONCURRENT (j = nx_min_2-1:nx_part_2, k = ny_min_2-1:ny_part_2, l = nz_min_2-1:nz_part_2)
	DO i = imin2, imax2
		uL2 (i,x_dir,j,k,l) = primL2(i,x_dir,j,k,l)*primL2(irho2,x_dir,j,k,l)
		uR2 (i,x_dir,j,k,l) = primR2(i,x_dir,j,k,l)*primR2(irho2,x_dir,j,k,l)
		IF(n_dim > 1) THEN
			uL2 (i,y_dir,j,k,l) = primL2(i,y_dir,j,k,l)*primL2(irho2,y_dir,j,k,l)
			uR2 (i,y_dir,j,k,l) = primR2(i,y_dir,j,k,l)*primR2(irho2,y_dir,j,k,l)
		END IF
		IF(n_dim > 2) THEN
			uL2 (i,z_dir,j,k,l) = primL2(i,z_dir,j,k,l)*primL2(irho2,z_dir,j,k,l)
			uR2 (i,z_dir,j,k,l) = primR2(i,z_dir,j,k,l)*primR2(irho2,z_dir,j,k,l)
		END IF
	END DO
	uL2 (irho2,x_dir,j,k,l) = primL2(irho2,x_dir,j,k,l)
	uR2 (irho2,x_dir,j,k,l) = primR2(irho2,x_dir,j,k,l)
	IF(n_dim > 1) THEN
		uL2 (irho2,y_dir,j,k,l) = primL2(irho2,y_dir,j,k,l)
		uR2 (irho2,y_dir,j,k,l) = primR2(irho2,y_dir,j,k,l)
	END IF
	IF(n_dim > 2) THEN
		uL2 (irho2,z_dir,j,k,l) = primL2(irho2,z_dir,j,k,l)
		uR2 (irho2,z_dir,j,k,l) = primR2(irho2,z_dir,j,k,l)
	END IF

	! NM energy equation !
	uL2 (itau2,x_dir,j,k,l) = primL2(irho2,x_dir,j,k,l)*(0.5D0*(primL2(ivel2_x,x_dir,j,k,l)**2 &
						    + primL2(ivel2_y,x_dir,j,k,l)**2 + primL2(ivel2_z,x_dir,j,k,l)**2) + eps2L(x_dir,j,k,l))
	uR2 (itau2,x_dir,j,k,l) = primR2(irho2,x_dir,j,k,l)*(0.5D0*(primR2(ivel2_x,x_dir,j,k,l)**2 &
							+ primR2(ivel2_y,x_dir,j,k,l)**2 + primR2(ivel2_z,x_dir,j,k,l)**2) + eps2R(x_dir,j,k,l))
	IF(n_dim > 1) THEN
		uL2 (itau2,y_dir,j,k,l) = primL2(irho2,y_dir,j,k,l)*(0.5D0*(primL2(ivel2_x,y_dir,j,k,l)**2 &
								+ primL2(ivel2_y,y_dir,j,k,l)**2 + primL2(ivel2_z,y_dir,j,k,l)**2) + eps2L(y_dir,j,k,l))
		uR2 (itau2,y_dir,j,k,l) = primR2(irho2,y_dir,j,k,l)*(0.5D0*(primR2(ivel2_x,y_dir,j,k,l)**2 &
								+ primR2(ivel2_y,y_dir,j,k,l)**2 + primR2(ivel2_z,y_dir,j,k,l)**2) + eps2R(y_dir,j,k,l))
	END IF
	IF(n_dim > 2) THEN
		uL2 (itau2,z_dir,j,k,l) = primL2(irho2,z_dir,j,k,l)*(0.5D0*(primL2(ivel2_x,z_dir,j,k,l)**2 &
								+ primL2(ivel2_y,z_dir,j,k,l)**2 + primL2(ivel2_z,z_dir,j,k,l)**2) + eps2L(z_dir,j,k,l))
		uR2 (itau2,z_dir,j,k,l) = primR2(irho2,z_dir,j,k,l)*(0.5D0*(primR2(ivel2_x,z_dir,j,k,l)**2 &
								+ primR2(ivel2_y,z_dir,j,k,l)**2 + primR2(ivel2_z,z_dir,j,k,l)**2) + eps2R(z_dir,j,k,l))
	END IF

	! dual energy !
	IF(dual_energy) THEN 
		uL2 (ieps2,x_dir,j,k,l) = primL2(ieps2,x_dir,j,k,l)
		uR2 (ieps2,x_dir,j,k,l) = primR2(ieps2,x_dir,j,k,l)
		IF(n_dim > 1) THEN
			uL2 (ieps2,y_dir,j,k,l) = primL2(ieps2,y_dir,j,k,l)
			uR2 (ieps2,y_dir,j,k,l) = primR2(ieps2,y_dir,j,k,l)
		END IF
		IF(n_dim > 2) THEN
			uL2 (ieps2,z_dir,j,k,l) = primL2(ieps2,z_dir,j,k,l)
			uR2 (ieps2,z_dir,j,k,l) = primR2(ieps2,z_dir,j,k,l)
		END IF
	END IF

	! For NM flux !
	DO i = imin2, imax2, 1
		fluxL2 (i,x_dir,j,k,l) = uL2 (i,x_dir,j,k,l) * primL2(ivel2_x,x_dir,j,k,l)
		fluxR2 (i,x_dir,j,k,l) = uR2 (i,x_dir,j,k,l) * primR2(ivel2_x,x_dir,j,k,l)
		IF(n_dim > 1) THEN
			fluxL2 (i,y_dir,j,k,l) = uL2 (i,y_dir,j,k,l) * primL2(ivel2_y,y_dir,j,k,l)
			fluxR2 (i,y_dir,j,k,l) = uR2 (i,y_dir,j,k,l) * primR2(ivel2_y,y_dir,j,k,l)
		END IF
		IF(n_dim > 2) THEN
			fluxL2 (i,z_dir,j,k,l) = uL2 (i,z_dir,j,k,l) * primL2(ivel2_z,z_dir,j,k,l)
			fluxR2 (i,z_dir,j,k,l) = uR2 (i,z_dir,j,k,l) * primR2(ivel2_z,z_dir,j,k,l)
		END IF
	ENDDO

	! Add the pressure term to x-momentum equation !
	fluxL2 (ivel2_x,x_dir,j,k,l) = fluxL2 (ivel2_x,x_dir,j,k,l) + primL2(itau2,x_dir,j,k,l)
	fluxR2 (ivel2_x,x_dir,j,k,l) = fluxR2 (ivel2_x,x_dir,j,k,l) + primR2(itau2,x_dir,j,k,l)
	IF(n_dim > 1) THEN
		fluxL2 (ivel2_y,y_dir,j,k,l) = fluxL2 (ivel2_y,y_dir,j,k,l) + primL2(itau2,y_dir,j,k,l)
		fluxR2 (ivel2_y,y_dir,j,k,l) = fluxR2 (ivel2_y,y_dir,j,k,l) + primR2(itau2,y_dir,j,k,l)
	END IF
	IF(n_dim > 2) THEN
		fluxL2 (ivel2_z,z_dir,j,k,l) = fluxL2 (ivel2_z,z_dir,j,k,l) + primL2(itau2,z_dir,j,k,l)
		fluxR2 (ivel2_z,z_dir,j,k,l) = fluxR2 (ivel2_z,z_dir,j,k,l) + primR2(itau2,z_dir,j,k,l)
	END IF

	! Add the presusre work done term to the energy equation             
	fluxL2 (itau2,x_dir,j,k,l) = fluxL2 (itau2,x_dir,j,k,l) + primL2(itau2,x_dir,j,k,l) * primL2(ivel2_x,x_dir,j,k,l)
	fluxR2 (itau2,x_dir,j,k,l) = fluxR2 (itau2,x_dir,j,k,l) + primR2(itau2,x_dir,j,k,l) * primR2(ivel2_x,x_dir,j,k,l) 
	IF(n_dim > 1) THEN
		fluxL2 (itau2,y_dir,j,k,l) = fluxL2 (itau2,y_dir,j,k,l) + primL2(itau2,y_dir,j,k,l) * primL2(ivel2_y,y_dir,j,k,l)
		fluxR2 (itau2,y_dir,j,k,l) = fluxR2 (itau2,y_dir,j,k,l) + primR2(itau2,y_dir,j,k,l) * primR2(ivel2_y,y_dir,j,k,l) 
	END IF
	IF(n_dim > 2) THEN
		fluxL2 (itau2,z_dir,j,k,l) = fluxL2 (itau2,z_dir,j,k,l) + primL2(itau2,z_dir,j,k,l) * primL2(ivel2_z,z_dir,j,k,l)
		fluxR2 (itau2,z_dir,j,k,l) = fluxR2 (itau2,z_dir,j,k,l) + primR2(itau2,z_dir,j,k,l) * primR2(ivel2_z,z_dir,j,k,l) 
	END IF

	! dual energy !
	IF(dual_energy) THEN  
		fluxL2 (ieps2,x_dir,j,k,l) = fluxL2 (ieps2,x_dir,j,k,l) + primL2(itau2,x_dir,j,k,l) * primL2(ivel2_x,x_dir,j,k,l)
		fluxR2 (ieps2,x_dir,j,k,l) = fluxR2 (ieps2,x_dir,j,k,l) + primR2(itau2,x_dir,j,k,l) * primR2(ivel2_x,x_dir,j,k,l)
		IF(n_dim > 1) THEN
			fluxL2 (ieps2,y_dir,j,k,l) = fluxL2 (ieps2,y_dir,j,k,l) + primL2(itau2,y_dir,j,k,l) * primL2(ivel2_y,y_dir,j,k,l)
			fluxR2 (ieps2,y_dir,j,k,l) = fluxR2 (ieps2,y_dir,j,k,l) + primR2(itau2,y_dir,j,k,l) * primR2(ivel2_y,y_dir,j,k,l)
		END IF
		IF(n_dim > 2) THEN
			fluxL2 (ieps2,z_dir,j,k,l) = fluxL2 (ieps2,z_dir,j,k,l) + primL2(itau2,z_dir,j,k,l) * primL2(ivel2_z,z_dir,j,k,l)
			fluxR2 (ieps2,z_dir,j,k,l) = fluxR2 (ieps2,z_dir,j,k,l) + primR2(itau2,z_dir,j,k,l) * primR2(ivel2_z,z_dir,j,k,l)
		END IF
	END IF
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE