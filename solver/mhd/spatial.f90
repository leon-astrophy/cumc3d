!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine prepares the data for spatial discretization,
! due ask the WENO_module to do the reconstruction
! and then combines the results for one Runge-Kutta sub-step
! Prototype developed by Wong Ka Wing in 2010 (or before?)
! Extended by Leung Shing Chi in 2016 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SPATIAL
USE DEFINITION
USE MHD_MODULE
USE PPM_MODULE
USE TVD_MODULE
USE WENO_MODULE
USE RIEMANN_MODULE
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! dummy variables !
INTEGER :: i, j, k, l, p

! Limits of density to be considered
REAL*8 :: rho_min1, rho_min2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find source term !

! Predefined source term !
CALL GET_SOURCE

! Custom source term !
CALL CUSTOM_SOURCE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do the x-sweep !

! Reconstruction !
IF(tvdvl_flag) THEN
	CALL tvdvl_reconx
ELSEIF(tvdmc_flag) THEN
	CALL tvdmc_reconx
ELSEIF(ppm1984_flag) THEN
	CALL ppm1984_reconx
ELSEIF(ppm2014_flag) THEN
	CALL ppm2014_reconx
ELSEIF(wenoz_flag) THEN
	CALL weno_reconx
END if

! Build states !
CALL BUILDSTATES(x_dir)

! For NM !
IF(LF_flag) THEN
	CALL LFNM(x_dir)
ELSEIF(HLL_flag) THEN
	CALL HLLNM(x_dir)
ELSEIF(HLLC_flag) THEN
	CALL HLLCNM(x_dir)
END IF

! Get source terms, fluxes, and flux graidents for NM
CALL FLUX_DIFF(x_dir)

! Get EMF !
CALL MHD_FLUX(x_dir)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do the y-sweep !
IF(n_dim > 1) THEN

	! Reconstruction !
	IF(tvdvl_flag) THEN
		CALL tvdvl_recony
	ELSEIF(tvdmc_flag) THEN
		CALL tvdmc_recony
	ELSEIF(ppm1984_flag) THEN
		CALL ppm1984_recony
	ELSEIF(ppm2014_flag) THEN
		CALL ppm2014_recony
	ELSEIF(wenoz_flag) THEN
		CALL weno_recony
	END if

		! Build states !
	CALL BUILDSTATES(y_dir)

	! For NM !
	IF(LF_flag) THEN
		CALL LFNM(y_dir)
	ELSEIF(HLL_flag) THEN
		CALL HLLNM(y_dir)
	ELSEIF(HLLC_flag) THEN
		CALL HLLCNM(y_dir)
	END IF

  ! Get source terms, fluxes, and flux graidents for NM
  CALL FLUX_DIFF(y_dir)

  ! Get EMF !
	CALL MHD_FLUX(y_dir)

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do the z-sweep !

IF(n_dim > 2) THEN

	! Reconstruction !
	IF(tvdvl_flag) THEN
		CALL tvdvl_reconz
	ELSEIF(tvdmc_flag) THEN
		CALL tvdmc_reconz
	ELSEIF(ppm1984_flag) THEN
		CALL ppm1984_reconz
	ELSEIF(ppm2014_flag) THEN
		CALL ppm2014_reconz
	ELSEIF(wenoz_flag) THEN
		CALL weno_reconz
	END if

	! Build states !
	CALL BUILDSTATES(z_dir)

	! For NM !
	IF(LF_flag) THEN
		CALL LFNM(z_dir)
	ELSEIF(HLL_flag) THEN
		CALL HLLNM(z_dir)
	ELSEIF(HLLC_flag) THEN
		CALL HLLCNM(z_dir)
	END IF

  ! Get source terms, fluxes, and flux graidents for NM
  CALL FLUX_DIFF(z_dir)

  ! Get EMF !
	CALL MHD_FLUX(z_dir)

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update Rungekutta operator !

! Constrained transport !
CALL flux_ct

! Get du/dt !
CALL l_operator

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine built the states for left and right edges for the horizontal directions !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GET_SOURCE
USE DEFINITION 
USE RIEMANN_MODULE
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k, l

! Pressure gradients !
REAL*8 :: dpdx2, dpdy2, dpdz2, bsquare

! Threshold for atmosphere density
REAL*8 :: rho_min1, rho_min2, factor, diff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now do the normal matter

! Set up the threshold
rho_min2 = 1.1D0 * prim2_a(irho2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Geometric sources

!$OMP PARALLEL PRIVATE(dpdx2, dpdy2, dpdz2, bsquare, diff, factor)

! Initialize !
!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
DO l = nz_min_2, nz_part_2
	DO k = ny_min_2, ny_part_2
		DO j = nx_min_2, nx_part_2
			sc2(imin2:imax2,j,k,l) = 0.0D0
		END DO
	END DO
END DO
!$OMP END DO

! Choose coordinate system, add geometric source term !
IF(coordinate_flag == 1) THEN
	!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
	DO l = nz_min_2, nz_part_2
  	DO k = ny_min_2, ny_part_2
   	 	DO j = nx_min_2, nx_part_2
			  diff = prim2(irho2,j,k,l) - rho_min2
      	factor = MAX(SIGN(1.0D0, diff), 0.0D0)
				bsquare = dot_product(prim2(ibx:ibz,j,k,l),prim2(ibx:ibz,j,k,l))
				sc2(ivel2_x,j,k,l) = sc2(ivel2_x,j,k,l) + (factor*prim2(itau2,j,k,l) + factor*prim2(irho2,j,k,l)*prim2(ivel2_y,j,k,l)**2 & 
													 											+ 0.5D0*bsquare - prim2(iby,j,k,l)**2)/x2(j)
				sc2(ivel2_y,j,k,l) = sc2(ivel2_y,j,k,l) - (factor*prim2(irho2,j,k,l)*prim2(ivel2_x,j,k,l)*prim2(ivel2_y,j,k,l) & 
													 											- prim2(ibx,j,k,l)*prim2(iby,j,k,l))/x2(j)
			END DO
		END DO
	END DO
	!$OMP END DO
ELSEIF(coordinate_flag == 2) THEN
	!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
	DO l = nz_min_2, nz_part_2
  	DO k = ny_min_2, ny_part_2
   	 	DO j = nx_min_2, nx_part_2
			  diff = prim2(irho2,j,k,l) - rho_min2
      	factor = MAX(SIGN(1.0D0, diff), 0.0D0)
				bsquare = dot_product(prim2(ibx:ibz,j,k,l),prim2(ibx:ibz,j,k,l))
				sc2(ivel2_x,j,k,l) = sc2(ivel2_x,j,k,l) + (2.0D0*factor*prim2(itau2,j,k,l) + prim2(ibx,j,k,l)**2 &
																								+ factor*prim2(irho2,j,k,l)*(prim2(ivel2_y,j,k,l)**2 + prim2(ivel2_z,j,k,l)**2))/x2(j)
				sc2(ivel2_y,j,k,l) = sc2(ivel2_y,j,k,l) + (factor*prim2(itau2,j,k,l) + factor*prim2(irho2,j,k,l)*prim2(ivel2_z,j,k,l)**2 & 
																								+ 0.5D0*bsquare - prim2(ibz,j,k,l)**2)*(COS(y2(k))/SIN(y2(k)))/x2(j) & 
																								- (factor*prim2(irho2,j,k,l)*prim2(ivel2_x,j,k,l)*prim2(ivel2_y,j,k,l) & 
																								- prim2(ibx,j,k,l)*prim2(iby,j,k,l))/x2(j)
				sc2(ivel2_z,j,k,l) = sc2(ivel2_z,j,k,l) - (factor*prim2(irho2,j,k,l)*prim2(ivel2_x,j,k,l)*prim2(ivel2_z,j,k,l) & 
																								- prim2(ibx,j,k,l)*prim2(ibz,j,k,l))/x2(j) &
																								- (factor*prim2(irho2,j,k,l)*prim2(ivel2_y,j,k,l)*prim2(ivel2_z,j,k,l) & 
																								- prim2(iby,j,k,l)*prim2(ibz,j,k,l))*(COS(y2(k))/SIN(y2(k)))/x2(j)
			END DO
		END DO
	END DO
	!$OMP END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! dual energy extra source !

IF(dual_energy) THEN

	! Choose according to coordinate system !
	IF(coordinate_flag == 0) THEN
		!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
		DO l = nz_min_2, nz_part_2
  		DO k = ny_min_2, ny_part_2
   	 		DO j = nx_min_2, nx_part_2
					CALL FINDGRADP(j,k,l,x_dir,dpdx2)
					CALL FINDGRADP(j,k,l,y_dir,dpdy2)
					CALL FINDGRADP(j,k,l,z_dir,dpdz2)
					sc2(ieps2,j,k,l) = sc2(ieps2,j,k,l) & 
													 + prim2(ivel2_x,j,k,l)*dpdx2 + prim2(ivel2_y,j,k,l)*dpdy2 + prim2(ivel2_z,j,k,l)*dpdz2
				END DO
			END DO
		END DO
		!$OMP END DO
	ELSEIF(coordinate_flag == 1) THEN
		!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
		DO l = nz_min_2, nz_part_2
  		DO k = ny_min_2, ny_part_2
   	 		DO j = nx_min_2, nx_part_2
					CALL FINDGRADP(j,k,l,x_dir,dpdx2)
					CALL FINDGRADP(j,k,l,y_dir,dpdy2)
					CALL FINDGRADP(j,k,l,z_dir,dpdz2)
					sc2(ieps2,j,k,l) = sc2(ieps2,j,k,l) & 
													 + prim2(ivel2_x,j,k,l)*dpdx2 + prim2(ivel2_y,j,k,l)*dpdy2/x2(j) + prim2(ivel2_z,j,k,l)*dpdz2
				END DO
			END DO
		END DO
		!$OMP END DO
	ELSEIF(coordinate_flag == 2) THEN
		!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
		DO l = nz_min_2, nz_part_2
  		DO k = ny_min_2, ny_part_2
   	 		DO j = nx_min_2, nx_part_2
					CALL FINDGRADP(j,k,l,x_dir,dpdx2)
					CALL FINDGRADP(j,k,l,y_dir,dpdy2)
					CALL FINDGRADP(j,k,l,z_dir,dpdz2)
					sc2(ieps2,j,k,l) = sc2(ieps2,j,k,l) & 
													 + prim2(ivel2_x,j,k,l)*dpdx2 + prim2(ivel2_y,j,k,l)*dpdy2/x2(j) + prim2(ivel2_z,j,k,l)*dpdz2/x2(j)/SIN(y2(k))
				END DO
			END DO
		END DO
		!$OMP END DO
	END IF

END IF
!$OMP END PARALLEL

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine built the states for left and right edges for the horizontal directions !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BUILDSTATES(dir_in)
USE DEFINITION 
USE MHD_MODULE
USE RIEMANN_MODULE
IMPLICIT NONE

! Integer !
INTEGER, INTENT(IN) :: dir_in

! Integer !
INTEGER :: ivn, ivt1, ivt2
INTEGER :: ibn, ibt1, ibt2
INTEGER :: kx, ky, kz
INTEGER :: i, j, k, l

! Dummy !
REAL*8 :: v2L, v2R, b2L, b2R, vbL, vbR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign !
kx = 0
ky = 0
kz = 0

! Assign !
IF(dir_in == x_dir) THEN
	ivn = ivel2_x
  ivt1 = ivel2_y
  ivt2 = ivel2_z
	ibn = ibx
	ibt1 = iby
	ibt2 = ibz
	kx = 1
ELSEIF(dir_in == y_dir) THEN
	ivn = ivel2_y
  ivt1 = ivel2_z
  ivt2 = ivel2_x
	ibn = iby
	ibt1 = ibz
	ibt2 = ibx
	ky = 1
ELSEIF(dir_in == z_dir) THEN
	ivn = ivel2_z
  ivt1 = ivel2_x
  ivt2 = ivel2_y
	ibn = ibz
	ibt1 = ibx
	ibt2 = iby
	kz = 1
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP PARALLEL PRIVATE(v2L, v2R, b2L, b2R, vbL, vbR)
!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
DO l = nz_min_2 - kz, nz_part_2
	DO k = ny_min_2 - ky, ny_part_2
		DO j = nx_min_2 - kx, nx_part_2

				! get dot product !
				b2L = dot_product(primL2(ibx:ibz,j,k,l), primL2(ibx:ibz,j,k,l))
				b2R = dot_product(primR2(ibx:ibz,j,k,l), primR2(ibx:ibz,j,k,l))
				v2L = dot_product(primL2(ivel2_x:ivel2_z,j,k,l), primL2(ivel2_x:ivel2_z,j,k,l))
				v2R = dot_product(primR2(ivel2_x:ivel2_z,j,k,l), primR2(ivel2_x:ivel2_z,j,k,l))
				vbL = dot_product(primL2(ibx:ibz,j,k,l), primL2(ivel2_x:ivel2_z,j,k,l))
				vbR = dot_product(primR2(ibx:ibz,j,k,l), primR2(ivel2_x:ivel2_z,j,k,l))

				! conservative variables !
				uL2 (imin2:imax2,j,k,l) = primL2(imin2:imax2,j,k,l)*primL2(irho2,j,k,l)
				uR2 (imin2:imax2,j,k,l) = primR2(imin2:imax2,j,k,l)*primR2(irho2,j,k,l)
				uL2 (irho2,j,k,l) = primL2(irho2,j,k,l)
				uR2 (irho2,j,k,l) = primR2(irho2,j,k,l)

				! Magnetic fields !
				uL2 (ibx:ibz,j,k,l) = primL2(ibx:ibz,j,k,l)
				uR2 (ibx:ibz,j,k,l) = primR2(ibx:ibz,j,k,l)

				! NM energy equation !
				uL2 (itau2,j,k,l) = primL2(irho2,j,k,l)*(0.5D0*v2L + eps2L(j,k,l)) + 0.5D0*b2L
				uR2 (itau2,j,k,l) = primR2(irho2,j,k,l)*(0.5D0*v2R + eps2R(j,k,l)) + 0.5D0*b2R

				! For NM flux !
				fluxL2 (imin2:imax2,j,k,l) = uL2 (imin2:imax2,j,k,l) * primL2(ivn,j,k,l)
				fluxR2 (imin2:imax2,j,k,l) = uR2 (imin2:imax2,j,k,l) * primR2(ivn,j,k,l)

				! Add the pressure term to x-momentum equation !
				fluxL2 (ivn,j,k,l) = fluxL2 (ivn,j,k,l) + primL2(itau2,j,k,l) + 0.5D0*(b2L) - primL2(ibn,j,k,l)**2
				fluxR2 (ivn,j,k,l) = fluxR2 (ivn,j,k,l) + primR2(itau2,j,k,l) + 0.5D0*(b2R) - primR2(ibn,j,k,l)**2

				! Other momentum fluxes !
				fluxL2(ivt1,j,k,l) = fluxL2(ivt1,j,k,l) - primL2(ibn,j,k,l)*primL2(ibt1,j,k,l)
				fluxR2(ivt1,j,k,l) = fluxR2(ivt1,j,k,l) - primR2(ibn,j,k,l)*primR2(ibt1,j,k,l)
				fluxL2(ivt2,j,k,l) = fluxL2(ivt2,j,k,l) - primL2(ibn,j,k,l)*primL2(ibt2,j,k,l)
				fluxR2(ivt2,j,k,l) = fluxR2(ivt2,j,k,l) - primR2(ibn,j,k,l)*primR2(ibt2,j,k,l)

				! Add the presusre work done term to the energy equation             
				fluxL2 (itau2,j,k,l) = fluxL2 (itau2,j,k,l) + (primL2(itau2,j,k,l) + 0.5D0*(b2L)) * primL2(ivn,j,k,l) - primL2(ibn,j,k,l)*vbL
				fluxR2 (itau2,j,k,l) = fluxR2 (itau2,j,k,l) + (primR2(itau2,j,k,l) + 0.5D0*(b2R)) * primR2(ivn,j,k,l) - primR2(ibn,j,k,l)*vbR

				! magentic field fluxes !
				fluxL2(ibn,j,k,l) = 0.0D0
				fluxR2(ibn,j,k,l) = 0.0D0
			  fluxL2(ibt1,j,k,l) = (primL2(ivn,j,k,l)*primL2(ibt1,j,k,l) - primL2(ivt1,j,k,l)*primL2(ibn,j,k,l))
				fluxR2(ibt1,j,k,l) = (primR2(ivn,j,k,l)*primR2(ibt1,j,k,l) - primR2(ivt1,j,k,l)*primR2(ibn,j,k,l))
			  fluxL2(ibt2,j,k,l) = (primL2(ivn,j,k,l)*primL2(ibt2,j,k,l) - primL2(ivt2,j,k,l)*primL2(ibn,j,k,l))
				fluxR2(ibt2,j,k,l) = (primR2(ivn,j,k,l)*primR2(ibt2,j,k,l) - primR2(ivt2,j,k,l)*primR2(ibn,j,k,l))

    END DO
  END DO
END DO
!$OMP END DO

! dual energy !
IF(dual_energy) THEN 
	!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
	DO l = nz_min_2 - kz, nz_part_2
		DO k = ny_min_2 - ky, ny_part_2
			DO j = nx_min_2 - kx, nx_part_2
				uL2 (ieps2,j,k,l) = primL2(ieps2,j,k,l)
				uR2 (ieps2,j,k,l) = primR2(ieps2,j,k,l)
				fluxL2 (ieps2,j,k,l) = uL2 (ieps2,j,k,l) * primL2(ivn,j,k,l)
				fluxR2 (ieps2,j,k,l) = uR2 (ieps2,j,k,l) * primR2(ivn,j,k,l)
				fluxL2 (ieps2,j,k,l) = fluxL2 (ieps2,j,k,l) + primL2(itau2,j,k,l) * primL2(ivn,j,k,l)
				fluxR2 (ieps2,j,k,l) = fluxR2 (ieps2,j,k,l) + primR2(itau2,j,k,l) * primR2(ivn,j,k,l)
    	END DO
  	END DO
	END DO
	!$OMP END DO
END IF

!$OMP END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine built the states for left and right edges for the horizontal directions !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FLUX_DIFF(dir_in)
USE DEFINITION 
USE RIEMANN_MODULE
IMPLICIT NONE

! Integer !
INTEGER, INTENT(IN) :: dir_in

! Integer !
INTEGER :: i, j, k, l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now calculate flux gradient, choose according to the coordinate system !
!$OMP PARALLEL
! x-direciton !
IF(dir_in == x_dir) THEN
	IF(coordinate_flag == 0) THEN
		!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
		DO l = nz_min_2, nz_part_2
			DO k = ny_min_2, ny_part_2
				DO j = nx_min_2, nx_part_2
					DO i = imin2, ibx - 1
						dflux_2 (i,j,k,l) = (flux_2 (i,j,k,l) - flux_2 (i,j-1,k,l)) / (dx2(j))
						l2(i,j,k,l) = - dflux_2 (i,j,k,l)
					END DO
				END DO
			END DO
		END DO
		!$OMP END DO
	ELSEIF(coordinate_flag == 1) THEN
		!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
		DO l = nz_min_2, nz_part_2
			DO k = ny_min_2, ny_part_2
				DO j = nx_min_2, nx_part_2
					DO i = imin2, ibx - 1
						dflux_2 (i,j,k,l) = (xF2(j)*flux_2 (i,j,k,l) - xF2(j-1)*flux_2 (i,j-1,k,l)) & 
															/ (x2(j)*dx2(j))
						l2(i,j,k,l) = - dflux_2 (i,j,k,l)
					END DO
				END DO
			END DO
		END DO
		!$OMP END DO
	ELSEIF(coordinate_flag == 2) THEN
		!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
		DO l = nz_min_2, nz_part_2
			DO k = ny_min_2, ny_part_2
				DO j = nx_min_2, nx_part_2
					DO i = imin2, ibx - 1
						dflux_2 (i,j,k,l) = (xF2(j)**2*flux_2 (i,j,k,l) - xF2(j-1)**2*flux_2 (i,j-1,k,l)) & 
															/ (x2(j)**2*dx2(j))
						l2(i,j,k,l) = - dflux_2 (i,j,k,l)
					END DO
				END DO
			END DO
		END DO
		!$OMP END DO
	END IF

! y-direciton !
ELSEIF(dir_in == y_dir) THEN
	IF(coordinate_flag == 0) THEN
		!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
		DO l = nz_min_2, nz_part_2
			DO k = ny_min_2, ny_part_2
				DO j = nx_min_2, nx_part_2
					DO i = imin2, ibx - 1
						dflux_2 (i,j,k,l) = (flux_2 (i,j,k,l) - flux_2 (i,j,k-1,l)) / (dy2(k))
						l2(i,j,k,l) = l2(i,j,k,l) - dflux_2 (i,j,k,l)
					END DO
				END DO
			END DO
		END DO
		!$OMP END DO
	ELSEIF(coordinate_flag == 1) THEN
		!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
		DO l = nz_min_2, nz_part_2
			DO k = ny_min_2, ny_part_2
				DO j = nx_min_2, nx_part_2
					DO i = imin2, ibx - 1
						dflux_2 (i,j,k,l) = (flux_2 (i,j,k,l) - flux_2 (i,j,k-1,l)) / (x2(j)*dy2(k))
						l2(i,j,k,l) = l2(i,j,k,l) - dflux_2 (i,j,k,l)
					END DO
				END DO
			END DO
		END DO
		!$OMP END DO
	ELSEIF(coordinate_flag == 2) THEN
		!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
		DO l = nz_min_2, nz_part_2
			DO k = ny_min_2, ny_part_2
				DO j = nx_min_2, nx_part_2
					DO i = imin2, ibx - 1
						dflux_2 (i,j,k,l) = (SIN(yF2(k))*flux_2 (i,j,k,l) - SIN(yF2(k-1))*flux_2 (i,j,k-1,l)) & 
															/ (x2(j)*SIN(y2(k))*dy2(k))
						l2(i,j,k,l) = l2(i,j,k,l) - dflux_2 (i,j,k,l)
					END DO
				END DO
			END DO
		END DO
		!$OMP END DO
	END IF

! z-direciton !
ELSEIF(dir_in == z_dir) THEN
	IF(coordinate_flag == 0) THEN
		!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
		DO l = nz_min_2, nz_part_2
			DO k = ny_min_2, ny_part_2
				DO j = nx_min_2, nx_part_2
					DO i = imin2, ibx - 1
						dflux_2 (i,j,k,l) = (flux_2 (i,j,k,l) - flux_2 (i,j,k,l-1)) / dz2(l)
						l2(i,j,k,l) = l2(i,j,k,l) - dflux_2 (i,j,k,l)
					END DO
				END DO
			END DO
		END DO
		!$OMP END DO
	ELSEIF(coordinate_flag == 1) THEN
		!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
		DO l = nz_min_2, nz_part_2
			DO k = ny_min_2, ny_part_2
				DO j = nx_min_2, nx_part_2
					DO i = imin2, ibx - 1
						dflux_2 (i,j,k,l) = (flux_2 (i,j,k,l) - flux_2 (i,j,k,l-1)) / dz2(l)
						l2(i,j,k,l) = l2(i,j,k,l) - dflux_2 (i,j,k,l)
					END DO
				END DO
			END DO
		END DO
		!$OMP END DO
	ELSEIF(coordinate_flag == 2) THEN
		!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
		DO l = nz_min_2, nz_part_2
			DO k = ny_min_2, ny_part_2
				DO j = nx_min_2, nx_part_2
					DO i = imin2, ibx - 1
						dflux_2 (i,j,k,l) = (flux_2 (i,j,k,l) - flux_2 (i,j,k,l-1)) & 
															/ (x2(j)*SIN(y2(k))*dz2(l))
						l2(i,j,k,l) = l2(i,j,k,l) - dflux_2 (i,j,k,l)
					END DO
				END DO
			END DO
		END DO
		!$OMP END DO
	END IF
END IF
!$OMP END PARALLEL 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine built the states for left and right edges for the horizontal directions !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE l_operator
USE DEFINITION 
USE MHD_MODULE
USE RIEMANN_MODULE
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k, l

! Real !
REAL*8 :: dbdt_l, dbdt_r

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP PARALLEL PRIVATE(dbdt_l, dbdt_r)
! Final step, get rungekutta operator, LHS of the hydro equation !
!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
DO l = nz_min_2, nz_part_2
	DO k = ny_min_2, ny_part_2
		DO j = nx_min_2, nx_part_2
			DO i = imin2, ibx - 1
				l2(i,j,k,l) = l2(i,j,k,l) + sc2(i,j,k,l)
			END DO
		END DO
	END DO
END DO
!$OMP END DO

! Magnetic field, special treatment for cylindrical coordinates !
IF(coordinate_flag == 0) THEN
	!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
	DO l = nz_min_2, nz_part_2
		DO k = ny_min_2, ny_part_2
			DO j = nx_min_2, nx_part_2

				! dbx/dt !
				dbdt_r = - ((efield_z(j,k,l) - efield_z(j,k-1,l))/(dy2(k)) - (efield_y(j,k,l) - efield_y(j,k,l-1))/(dz2(l)))
				dbdt_l = - ((efield_z(j-1,k,l) - efield_z(j-1,k-1,l))/(dy2(k)) - (efield_y(j-1,k,l) - efield_y(j-1,k,l-1))/(dz2(l)))
				lint_x(j,k,l) = dbdt_r
				lint_x(j-1,k,l) = dbdt_l
				l2(ibx,j,k,l) = 0.5D0*(dbdt_l + dbdt_r)

				! dby/dt !
				dbdt_r = - ((efield_x(j,k,l) - efield_x(j,k,l-1))/(dz2(l)) - (efield_z(j,k,l) - efield_z(j-1,k,l))/(dx2(j)))
				dbdt_l = - ((efield_x(j,k-1,l) - efield_x(j,k-1,l-1))/(dz2(l)) - (efield_z(j,k-1,l) - efield_z(j-1,k-1,l))/(dx2(j)))
				lint_y(j,k,l) = dbdt_r
				lint_y(j,k-1,l) = dbdt_l
				l2(iby,j,k,l) = 0.5D0*(dbdt_l + dbdt_r)

				! dbz/dt !
				dbdt_r = - ((efield_y(j,k,l) - efield_y(j-1,k,l))/(dx2(j)) - (efield_x(j,k,l) - efield_x(j,k-1,l))/(dy2(k)))
				dbdt_l = - ((efield_y(j,k,l-1) - efield_y(j-1,k,l-1))/(dx2(j)) - (efield_x(j,k,l-1) - efield_x(j,k-1,l-1))/(dy2(k)))
				lint_z(j,k,l) = dbdt_r
				lint_z(j,k,l-1) = dbdt_l
				l2(ibz,j,k,l) = 0.5D0*(dbdt_l + dbdt_r)

			END DO
		END DO
	END DO
	!$OMP END DO
ELSEIF(coordinate_flag == 1) THEN
	!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
	DO l = nz_min_2, nz_part_2
		DO k = ny_min_2, ny_part_2
			DO j = nx_min_2, nx_part_2

				! dbx/dt !
				dbdt_r = - ((efield_z(j,k,l) - efield_z(j,k-1,l))/(xF2(j)*dy2(k)) - (efield_y(j,k,l) - efield_y(j,k,l-1))/(dz2(l)))
				dbdt_l = - ((efield_z(j-1,k,l) - efield_z(j-1,k-1,l))/(xF2(j-1)*dy2(k)) - (efield_y(j-1,k,l) - efield_y(j-1,k,l-1))/(dz2(l)))
				lint_x(j,k,l) = dbdt_r
				lint_x(j-1,k,l) = dbdt_l
				l2(ibx,j,k,l) = 0.5D0*(dbdt_l + dbdt_r)

				! dby/dt !
				dbdt_r = - ((efield_x(j,k,l) - efield_x(j,k,l-1))/(dz2(l)) - (efield_z(j,k,l) - efield_z(j-1,k,l))/(dx2(j)))
				dbdt_l = - ((efield_x(j,k-1,l) - efield_x(j,k-1,l-1))/(dz2(l)) - (efield_z(j,k-1,l) - efield_z(j-1,k-1,l))/(dx2(j)))
				lint_y(j,k,l) = dbdt_r
				lint_y(j,k-1,l) = dbdt_l
				l2(iby,j,k,l) = 0.5D0*(dbdt_l + dbdt_r)

				! dbz/dt !
				dbdt_r = - ((xF2(j)*efield_y(j,k,l) - xF2(j-1)*efield_y(j-1,k,l))/(dx2(j)) - (efield_x(j,k,l) - efield_x(j,k-1,l))/(dy2(k)))/x2(j)
				dbdt_l = - ((xF2(j)*efield_y(j,k,l-1) - xF2(j-1)*efield_y(j-1,k,l-1))/(dx2(j)) - (efield_x(j,k,l-1) - efield_x(j,k-1,l-1))/(dy2(k)))/x2(j)
				lint_z(j,k,l) = dbdt_r
				lint_z(j,k,l-1) = dbdt_l
				l2(ibz,j,k,l) = 0.5D0*(dbdt_l + dbdt_r)

			END DO
		END DO
	END DO
	!$OMP END DO
ELSEIF(coordinate_flag == 2) THEN
	!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
	DO l = nz_min_2, nz_part_2
		DO k = ny_min_2, ny_part_2
			DO j = nx_min_2, nx_part_2

				! dbx/dt !
				dbdt_r = - ((SIN(yF2(k))*efield_z(j,k,l) - SIN(yF2(k-1))*efield_z(j,k-1,l))/(dy2(k)) - (efield_y(j,k,l) - efield_y(j,k,l-1))/(dz2(l)))/(xF2(j)*SIN(y2(k)))
				dbdt_l = - ((SIN(yF2(k))*efield_z(j-1,k,l) - SIN(yF2(k-1))*efield_z(j-1,k-1,l))/(dy2(k)) - (efield_y(j-1,k,l) - efield_y(j-1,k,l-1))/(dz2(l)))/(xF2(j-1)*SIN(y2(k)))
				lint_x(j,k,l) = dbdt_r
				lint_x(j-1,k,l) = dbdt_l
				l2(ibx,j,k,l) = 0.5D0*(dbdt_l + dbdt_r)

				! dby/dt !
				dbdt_r = - ((efield_x(j,k,l) - efield_x(j,k,l-1))/(SIN(yF2(k))*dz2(l)) - (xF2(j)*efield_z(j,k,l) - xF2(j-1)*efield_z(j-1,k,l))/(dx2(j)))/x2(j)
				dbdt_l = - ((efield_x(j,k-1,l) - efield_x(j,k-1,l-1))/(SIN(yF2(k-1))*dz2(l)) - (xF2(j)*efield_z(j,k-1,l) - xF2(j-1)*efield_z(j-1,k-1,l))/(dx2(j)))/x2(j)
				lint_y(j,k,l) = dbdt_r
				lint_y(j,k-1,l) = dbdt_l
				l2(iby,j,k,l) = 0.5D0*(dbdt_l + dbdt_r)

				! dbz/dt !
				dbdt_r = - ((xF2(j)*efield_y(j,k,l) - xF2(j-1)*efield_y(j-1,k,l))/(dx2(j)) - (efield_x(j,k,l) - efield_x(j,k-1,l))/(dy2(k)))/x2(j)
				dbdt_l = - ((xF2(j)*efield_y(j,k,l-1) - xF2(j-1)*efield_y(j-1,k,l-1))/(dx2(j)) - (efield_x(j,k,l-1) - efield_x(j,k-1,l-1))/(dy2(k)))/x2(j)
				lint_z(j,k,l) = dbdt_r
				lint_z(j,k,l-1) = dbdt_l
				l2(ibz,j,k,l) = 0.5D0*(dbdt_l + dbdt_r)
			
			END DO
		END DO
	END DO
	!$OMP END DO
END IF
!$OMP END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE