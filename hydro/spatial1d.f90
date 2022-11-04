!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine prepares the data for spatial discretization,
! due ask the WENO_module to do the reconstruction
! and then combines the results for one Runge-Kutta sub-step
! Prototype developed by Wong Ka Wing in 2010 (or before?)
! Extended by Leung Shing Chi in 2016 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SPATIAL1D
USE DEFINITION
IMPLICIT NONE

! dummy variables !
INTEGER :: i, j, k

! Limits of density to be considered
REAL (DP) :: rho_min1, rho_min2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Find frame velocity 
IF (movinggriddm_flag == 1) THEN
	!CALL FINDFRAMEVEL_DM
END IF
IF (movinggridnm_flag == 1) THEN
	!CALL FINDFRAMEVEL_NM
END IF

! Initialize source term and fluxes !
IF(RUNDM_flag == 1) THEN
	src1 = 0.0D0
	flx_1 = 0.0D0
	dfdx1 = 0.0D0
END IF

! For NM !
src2 = 0.0D0
flx_2 = 0.0D0
dfdx2 = 0.0D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Find effective speed !
CALL FINDALPHA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section 1: Reconstructions and build states in the x-directions !

! Reconstruct using WENO !
CALL WENO_X

! Pressure graident !
IF(dual_energy == 1) THEN
	CALL FINDDPDX
END IF

! Build the fluxes and states !
CALL BUILDSTATES_X

! Choose an appropriate riemann solver !
IF(RUNDM_flag == 1) THEN
	CALL LFDM_X(flx_1)
END IF

! For NM !
CALL LFNM_X(flx_2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section 2: Get source terms

! For DM sector
IF(runDM_flag == 1) THEN
	rho_min1 = 1.1D0 * prim1_a(irho1)

	! Gravitational source term !
	DO j = nx_min_1, nx_part_1
		DO k = ny_min_1, ny_part_1
			DO l = nz_min_1, nz_part_1
				IF(prim1(j,k,l,irho1) > rho_min1) THEN
					IF(lapse_flag) THEN
						src1(j,k,l,ivel1_x) = alapse_1(j,k,l) * (prim1(j,k,l,irho1) - p1(j,k,l)) * phi1_x(j,k,l)
					ELSE
						src1(j,k,l,ivel1_x) = prim1(j,k,l,irho1) * phi1_x(j,k,l)
					END IF
				END IF
			END DO
		END DO
	ENDDO

	! Add source term to momentum equation according to the coordinate system !
	IF(coordinate_flag == 1) THEN
		sa1(:,:,ivel1_r) = - p1(:,:) - rho1(:,:)*vel1_p(:,:)**2
	ELSEIF(coordinate_flag == 2) THEN
	   	IF(lapse_flag == 1) THEN 
			sa1(:,:,ivel1_r) = - 2.0D0*alapse_1(:,:)*p1(:,:) - rho1(:,:)*alapse_1(:,:)*(vel1_z(:,:)**2 + vel1_p(:,:)**2)	
	   	ELSE
			sa1(:,:,ivel1_r) = - 2.0D0*p1(:,:) - rho1(:,:)*(vel1_z(:,:)**2 + vel1_p(:,:)**2)
	   END IF	
	ELSEIF(coordinate_flag == 3) THEN
	END IF

	! Moving grid !
	IF(movinggriddm_flag == 1) THEN
		sb1 (:,:,:) = sb1 (:,:,:) + u1 (:,:,:) * (3.0D0 * vel1_max / radius1)
	END IF

END IF

! Threshold density !
rho_min2 = 1.1D0 * prim2_a(irho1)

! For NM sector
	DO j = nx_min_1, nx_part_1
		DO k = ny_min_1, ny_part_1
			DO l = nz_min_1, nz_part_1
		IF(rho2(j,k) > rho_min2) THEN
			IF(lapse_flag == 1) THEN
				sb2(j,k,ivel2_r) = alapse_2(j,k) * (rho2(j,k) - p2(j,k)) * phi2_r(j,k)
				IF(coordinate_flag == 2) THEN
					sb2(j,k,itau2) = alapse_2(j,k)*rho2(j,k)*(vel2_r(j,k)*phi2_r(j,k) + vel2_z(j,k)*phi2_z(j,k)/r2(j))
				ELSE
				  	sb2(j,k,itau2) = alapse_2(j,k)*rho2(j,k)*(vel2_r(j,k)*phi2_r(j,k) + vel2_z(j,k)*phi2_z(j,k))
				END IF
	 	   	ELSE
				sb2(j,k,ivel2_r) = rho2(j,k) * phi2_r(j,k)
				IF(coordinate_flag == 2) THEN
				  	sb2(j,k,itau2) = rho2(j,k)*(vel2_r(j,k)*phi2_r(j,k) + vel2_z(j,k)*phi2_z(j,k)/r2(j))
				ELSE
				  	sb2(j,k,itau2) = rho2(j,k)*(vel2_r(j,k)*phi2_r(j,k) + vel2_z(j,k)*phi2_z(j,k))
				END IF
				END IF
		   END IF
		END IF
	END DO
END DO

IF(dual_energy == 1) THEN
   	IF(lapse_flag == 1) THEN
		IF(coordinate_flag == 2) THEN
			sb2(:,:,ieps2) = - (vel2_r(:,:)*dp2dr(:,:) + vel2_z(:,:)*dp2dz(:,:)/rad2(:,:)) + alapse_2(:,:)*p2(:,:)*(vel2_r(:,:)*phi2_r(:,:) + vel2_z(:,:)*phi2_z(:,:)/rad2(:,:))
		ELSE
			sb2(:,:,ieps2) = - (vel2_r(:,:)*dp2dr(:,:) + vel2_z(:,:)*dp2dz(:,:)) + alapse_2(:,:)*p2(:,:)*(vel2_r(:,:)*phi2_r(:,:) + vel2_z(:,:)*phi2_z(:,:))
		END IF
   	ELSE
		IF(coordinate_flag == 2) THEN
			sb2(:,:,ieps2) = - (vel2_r(:,:)*dp2dr(:,:) + vel2_z(:,:)*dp2dz(:,:)/rad2(:,:))
		ELSE
			sb2(:,:,ieps2) = - (vel2_r(:,:)*dp2dr(:,:) + vel2_z(:,:)*dp2dz(:,:))
		END IF
   END IF
END IF

! Choose coordinate system 
IF(coordinate_flag == 1) THEN
	sa2(:,:,ivel2_r) = - p2(:,:) - rho2(:,:)*vel2_p(:,:)**2
ELSEIF(coordinate_flag == 2) THEN
   	IF(lapse_flag == 1) THEN
		sa2(:,:,ivel2_r) = - 2.0D0*alapse_2(:,:)*p2(:,:) - rho2(:,:)*alapse_2(:,:)*(vel2_z(:,:)**2 + vel2_p(:,:)**2)
   	ELSE
		sa2(:,:,ivel2_r) = - 2.0D0*p2(:,:) - rho2(:,:)*(vel2_z(:,:)**2 + vel2_p(:,:)**2)
   	END IF
END IF

! Moving grid !
IF(movinggridnm_flag == 1) THEN
	sb2 (:,:,:) = sb2 (:,:,:) + u2 (:,:,:) * (3.0D0 * vel2_max / radius2)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section 4: Get flux derivatives

! Now calculate dfdx accordingly to the corrected flux, do it case by case
IF(coordinate_flag == 0) THEN
	IF(RUNDM_flag == 1) THEN
		DO j = 1, length_step_r_part_1, 1
			dfdr1 (j, :, :) = (fluxr_1 (j, :, :) - fluxr_1 (j - 1, :, :)) / dr1(j)
		END DO
	END IF

	! Do for NM !
	DO j = 1, length_step_r_part_2, 1
		dfdr2 (j, :, :) = (fluxr_2 (j, :, :) - fluxr_2 (j - 1, :, :)) / dr2(j)
	END DO

ELSEIF (coordinate_flag == 1) THEN
	IF(RUNDM_flag == 1) THEN
		DO j = 1, length_step_r_part_1, 1
			dfdr1 (j, :, :) = (rF1(j)*fluxr_1 (j, :, :) - rF1(j-1)*fluxr_1 (j - 1, :, :)) / (r1(j)*dr1(j))
		END DO

	END IF

	! Do for NM !
	DO j = 1, length_step_r_part_2, 1
		dfdr2 (j, :, :) = (rF2(j)*fluxr_2 (j, :, :) - rF2(j-1)*fluxr_2 (j - 1, :, :)) / (r2(j)*dr2(j))
	END DO

ELSEIF (coordinate_flag == 2) THEN
	IF(RUNDM_flag == 1) THEN
		DO j = 1, length_step_r_part_1, 1
			dfdr1 (j, :, :) = (rF1(j)**2*fluxr_1 (j, :, :) - rF1(j-1)**2*fluxr_1 (j - 1, :, :)) / (r1(j)**2*dr1(j))
		END DO

	END IF

	! Do for NM !
	DO j = 1, length_step_r_part_2, 1
		dfdr2 (j, :, :) = (rF2(j)**2*fluxr_2 (j, :, :) - rF2(j-1)**2*fluxr_2 (j - 1, :, :)) / (r2(j)**2*dr2(j))
	END DO

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section 5: Calculate l (SSPRK54)

! Initialization
l1 = 0.0D0
l2 = 0.0D0

! assign !
IF(coordinate_flag == 0) THEN
	IF(runDM_flag == 1) THEN
		l1(:,:,:) = - dfdx1 (:,:,:) - sb1 (:,:,:) 
   	ENDIF
   	l2(:,:,:) = - dfdr2 (:,:,:) - sb2 (:,:,:)
ELSEIF(coordinate_flag == 1 .OR. coordinate_flag == 2) THEN
   	IF(runDM_flag == 1) THEN
        DO j = 1, length_step_r_part_1, 1
            l1(j,:,:) = - dfdr1 (j,:,:) - sb1 (j,:,:) - sa1 (j,:,:) / r1(j)
        ENDDO
   	ENDIF
   	DO j = 1, length_step_r_part_2, 1
        l2(j,:,:) = - dfdr2 (j,:,:) - sb2 (j,:,:) - sa2 (j,:,:) / r2(j)
   	ENDDO
ENDIF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine built the states for left and right edges for the horizontal directions !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BUILDSTATES_X
USE DEFINITION 
USE RIEMANN_MODULE
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Build conserved variables !
If (RUNDM_flag == 1) then
	DO i = imin1, imax1
		uL1 (:, :, i) = primL1(:, :, i)*primL1(:, :, irho1)
		uR1 (:, :, i) = primR1(:, :, i)*primR1(:, :, irho1)
	END DO
	uL1 (:, :, irho1) = primL1(:, :, irho1)
	uR1 (:, :, irho1) = primR1(:, :, irho1)
END IF

! Do the same for NM !
DO i = imin2, imax2
	uL2 (:, :, i) = primL2(:, :, i)*primL2(:, :, irho2)
	uR2 (:, :, i) = primR2(:, :, i)*primR2(:, :, irho2)
END DO
uL2 (:, :, irho2) = primL2(:, :, irho2)
uR2 (:, :, irho2) = primR2(:, :, irho2)

! NM energy equation !
IF(nm_epsilon == 1) THEN
	uL2 (:, :, itau2) = primL2(:, :, irho2)*(0.5D0*(primL2(:, :, ivel2_r)**2 + primL2(:, :, ivel2_z)**2) + eps2L(:,:))
	uR2 (:, :, itau2) = primR2(:, :, irho2)*(0.5D0*(primR2(:, :, ivel2_r)**2 + primR2(:, :, ivel2_z)**2) + eps2R(:,:))
	IF(rotationnm_flag == 1) THEN
		uL2 (:, :, itau2) = uL2 (:, :, itau2) + primL2(:, :, irho2)*0.5D0*primL2(:, :, ivel2_p)**2
		uR2 (:, :, itau2) = uR2 (:, :, itau2) + primR2(:, :, irho2)*0.5D0*primR2(:, :, ivel2_p)**2
	END IF
END IF
IF(dual_energy == 1) THEN 
	uL2 (:, :, ieps2) = primL2(:, :, ieps2)
	uR2 (:, :, ieps2) = primR2(:, :, ieps2)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Build fluxes For DM !
If (RUNDM_flag == 1) then
	DO i = imin1, imax1
		fluxL1 (:, :, i) = uL1 (:, :, i) * primL1(:, :, ivel1_r)
		fluxR1 (:, :, i) = uR1 (:, :, i) * primR1(:, :, ivel1_r)
	END DO

   	! Add the pressure term to r-momentum equation !
    fluxL1 (:,:,ivel1_r) = fluxL1 (:,:,ivel1_r) + p1L(:,:)
	fluxR1 (:,:,ivel1_r) = fluxR1 (:,:,ivel1_r) + p1R(:,:)

	! Modify the flux !
	IF(lapse_flag == 1) THEN 
   	   	DO i = imin1, imax1, 1
			fluxL1 (:,:,i) = fluxL1 (:,:,i) * ap1L(:,:)
			fluxR1 (:,:,i) = fluxR1 (:,:,i) * ap1R(:,:)
	   	ENDDO
	END IF

	! Moving grid !
	If(movinggriddm_flag == 1) THEN
		DO i = imin1, imax1, 1
         	fluxL1 (:,:,i) = fluxL1 (:,:,i) - uL1(:,:,i)*vf1rL(:,:)
			fluxR1 (:,:,i) = fluxR1 (:,:,i) - uR1(:,:,i)*vf1rR(:,:)
		END DO
	END IF
END IF

! For NM !
DO i = imin2, imax2, 1
	fluxL2 (:, :, i) = uL2 (:, :, i) * primL2(:, :, ivel2_r)
	fluxR2 (:, :, i) = uR2 (:, :, i) * primR2(:, :, ivel2_r)
ENDDO

! Add the pressure term to r-momentum equation !
fluxL2 (:,:,ivel2_r) = fluxL2 (:,:,ivel2_r) + p2L(:,:)
fluxR2 (:,:,ivel2_r) = fluxR2 (:,:,ivel2_r) + p2R(:,:)

! Add the presusre work done term to the energy equation             
IF(nm_epsilon == 1) THEN  
	fluxL2 (:,:,itau2) = fluxL2 (:,:,itau2) + p2L(:,:) * primL2(:, :, ivel2_r)
	fluxR2 (:,:,itau2) = fluxR2 (:,:,itau2) + p2R(:,:) * primR2(:, :, ivel2_r)
END IF  
IF(dual_energy == 1) THEN  
	fluxL2 (:,:,ieps2) = fluxL2 (:,:,ieps2) + p2L(:,:) * primL2(:, :, ivel2_r)
	fluxR2 (:,:,ieps2) = fluxR2 (:,:,ieps2) + p2R(:,:) * primR2(:, :, ivel2_r)
END IF

! Modify !
IF(lapse_flag == 1) THEN
   DO i = imin2, imax2, 1
		fluxL2 (:,:,i) = fluxL2 (:,:,i) * ap2L (:,:)
		fluxR2 (:,:,i) = fluxR2 (:,:,i) * ap2R (:,:)
   ENDDO
END IF

! Extra flux term for moving grid !
If(movinggridnm_flag == 1) THEN
	DO i = imin2, imax2, 1
        fluxL2 (:,:,i) = fluxL2 (:,:,i) - uL2(:,:,i)*vf2rL(:,:)
		fluxR2 (:,:,i) = fluxR2 (:,:,i) - uR2(:,:,i)*vf2rR(:,:)
	END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE