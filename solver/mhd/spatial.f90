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
USE PPMC_MODULE
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
ELSEIF(ppm_flag) THEN
	CALL ppm_reconx
ELSEIF(ppmc_flag) THEN
	CALL ppmc_reconx
END IF

! Build states !
CALL BUILDSTATES(x_dir)

! For NM !
IF(LF_flag) THEN
	CALL LFNM(x_dir)
ELSEIF(HLL_flag) THEN
	CALL HLLNM(x_dir)
ELSEIF(HLLC_flag) THEN
	CALL HLLCNM(x_dir)
ELSEIF(HLLD_flag) THEN
	CALL HLLDNM(x_dir)
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
	ELSEIF(ppm_flag) THEN
		CALL ppm_recony
	ELSEIF(ppmc_flag) THEN
		CALL ppmc_recony
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
	ELSEIF(HLLD_flag) THEN
		CALL HLLDNM(y_dir)
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
	ELSEIF(ppm_flag) THEN
		CALL ppm_reconz
	ELSEIF(ppmc_flag) THEN
		CALL ppmc_reconz
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
	ELSEIF(HLLD_flag) THEN
		CALL HLLDNM(z_dir)
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

! Check timing with or without openmp
INTEGER :: time_start, time_end
INTEGER :: cr
REAL*8 :: rate
CALL system_clock(count_rate=cr)
rate = REAL(cr)

CALL system_clock(time_start)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now do the normal matter

! Geometric sources
!$OMP PARALLEL PRIVATE(dpdx2, dpdy2, dpdz2, bsquare, diff, factor)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize !
!$OMP DO COLLAPSE(4) SCHEDULE(STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT) 
DO l = nz_min_2, nz_part_2
	DO k = ny_min_2, ny_part_2
		DO j = nx_min_2, nx_part_2
			DO i = imin2, imax2
				sc2(i,j,k,l) = 0.0D0
			END DO
		END DO
	END DO
END DO
!$ACC END PARALLEL
!$OMP END DO

! Choose coordinate system, add geometric source term !
IF(coordinate_flag == 1) THEN
	!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
	!$ACC PARALLEL  LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) PRIVATE(bsquare, diff, factor)
	DO l = nz_min_2, nz_part_2
  	DO k = ny_min_2, ny_part_2
   	 	DO j = nx_min_2, nx_part_2
				rho_min2 = 1.1D0 * prim2_a(irho2)
			  diff = prim2(irho2,j,k,l) - rho_min2
      	factor = MAX(SIGN(1.0D0, diff), 0.0D0)
				bsquare = dot_product(prim2(ibx:ibz,j,k,l),prim2(ibx:ibz,j,k,l))
				sc2(ivel2_x,j,k,l) = sc2(ivel2_x,j,k,l) + (factor*prim2(itau2,j,k,l) + factor*prim2(irho2,j,k,l)*prim2(ivel2_y,j,k,l)**2 & 
													 											+ 0.5D0*bsquare - prim2(iby,j,k,l)**2)*(2.0D0*dx2(j)/dx2_sq(j))!/x2(j)	
																																											
				sc2(ivel2_y,j,k,l) = sc2(ivel2_y,j,k,l) - (factor*prim2(irho2,j,k,l)*prim2(ivel2_x,j,k,l)*prim2(ivel2_y,j,k,l) & 
													 											- prim2(ibx,j,k,l)*prim2(iby,j,k,l))*(2.0D0*dx2(j)/dx2_sq(j))!/x2(j)		
			END DO
		END DO
	END DO
	!$ACC END PARALLEL
	!$OMP END DO
ELSEIF(coordinate_flag == 2) THEN
	!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
	!$ACC PARALLEL  LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) PRIVATE(bsquare, diff, factor)
	DO l = nz_min_2, nz_part_2
  	DO k = ny_min_2, ny_part_2
   	 	DO j = nx_min_2, nx_part_2
				rho_min2 = 1.1D0 * prim2_a(irho2)
			  diff = prim2(irho2,j,k,l) - rho_min2
      	factor = MAX(SIGN(1.0D0, diff), 0.0D0)
				bsquare = dot_product(prim2(ibx:ibz,j,k,l),prim2(ibx:ibz,j,k,l))
				sc2(ivel2_x,j,k,l) = sc2(ivel2_x,j,k,l) + (2.0D0*factor*prim2(itau2,j,k,l) + prim2(ibx,j,k,l)**2 &
																								+ factor*prim2(irho2,j,k,l)*(prim2(ivel2_y,j,k,l)**2 + prim2(ivel2_z,j,k,l)**2)) & 
																								*(1.5D0*dx2_sq(j)/dx2_cb(j)) !/x2(j)

				sc2(ivel2_y,j,k,l) = sc2(ivel2_y,j,k,l) + (factor*prim2(itau2,j,k,l) + factor*prim2(irho2,j,k,l)*prim2(ivel2_z,j,k,l)**2 & 
																								+ 0.5D0*bsquare - prim2(ibz,j,k,l)**2)*(1.5D0*dx2_sq(j)/dx2_cb(j))*(dsin2(k)/dcos2(k)) & !/TAN(y2(k)) & !/x2(j)
																						
																								- (factor*prim2(irho2,j,k,l)*prim2(ivel2_x,j,k,l)*prim2(ivel2_y,j,k,l) & 
																								- prim2(ibx,j,k,l)*prim2(iby,j,k,l))*(1.5D0*dx2_sq(j)/dx2_cb(j)) !/x2(j) 

				sc2(ivel2_z,j,k,l) = sc2(ivel2_z,j,k,l) - (factor*prim2(irho2,j,k,l)*prim2(ivel2_x,j,k,l)*prim2(ivel2_z,j,k,l) & 
																								- prim2(ibx,j,k,l)*prim2(ibz,j,k,l))*(1.5D0*dx2_sq(j)/dx2_cb(j)) & !/x2(j)

																							  - (factor*prim2(irho2,j,k,l)*prim2(ivel2_y,j,k,l)*prim2(ivel2_z,j,k,l) & 
																								- prim2(iby,j,k,l)*prim2(ibz,j,k,l))*(1.5D0*dx2_sq(j)/dx2_cb(j))*(dsin2(k)/dcos2(k)) !/TAN(y2(k)) !/x2(j)
			END DO
		END DO
	END DO
	!$ACC END PARALLEL
	!$OMP END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! dual energy extra source !

! Choose according to coordinate system ! 
IF(dual_energy) THEN
	IF(coordinate_flag == 0) THEN
		!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
		!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) PRIVATE(dpdx2, dpdy2, dpdz2)
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
		!$ACC END PARALLEL
		!$OMP END DO
	ELSEIF(coordinate_flag == 1) THEN
		!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
		!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) PRIVATE(dpdx2, dpdy2, dpdz2)
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
		!$ACC END PARALLEL
		!$OMP END DO
	ELSEIF(coordinate_flag == 2) THEN
		!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
		!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) PRIVATE(dpdx2, dpdy2, dpdz2)
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
		!$ACC END PARALLEL
		!$OMP END DO
	END IF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL system_clock(time_end)
#ifdef DEBUG
WRITE(*,*) 'getsource = ', REAL(time_end - time_start) / rate
#endif

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the pressure gradient !                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDGRADP(j_in,k_in,l_in,dir_in,dp_out)
USE DEFINITION
USE RIEMANN_MODULE
IMPLICIT NONE

! Integer !
INTEGER, INTENT(IN) :: dir_in

! Integer !
INTEGER, INTENT(IN) :: j_in,k_in,l_in

! Real !
REAL*8, INTENT(OUT) :: dp_out

! Case by case !
IF(dir_in == x_dir) THEN
ELSEIF(dir_in == y_dir) THEN
ELSEIF(dir_in == z_dir) THEN
END IF

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

! Check timing with or without openmp
INTEGER :: time_start, time_end
INTEGER :: cr
REAL*8 :: rate

CALL system_clock(count_rate=cr)
rate = REAL(cr)

CALL system_clock(time_start)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign !
kx = 1
ky = 1
kz = 1

! Assign !
IF(dir_in == x_dir) THEN
	ivn = ivel2_x
  ivt1 = ivel2_y
  ivt2 = ivel2_z
	ibn = ibx
	ibt1 = iby
	ibt2 = ibz
	kx = 0
ELSEIF(dir_in == y_dir) THEN
	ivn = ivel2_y
  ivt1 = ivel2_z
  ivt2 = ivel2_x
	ibn = iby
	ibt1 = ibz
	ibt2 = ibx
	ky = 0
ELSEIF(dir_in == z_dir) THEN
	ivn = ivel2_z
  ivt1 = ivel2_x
  ivt2 = ivel2_y
	ibn = ibz
	ibt1 = ibx
	ibt2 = iby
	kz = 0
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP PARALLEL PRIVATE(v2L, v2R, b2L, b2R, vbL, vbR)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get flux !

!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE (3) DEFAULT(PRESENT) PRIVATE(v2L, v2R, b2L, b2R, vbL, vbR)
DO l = nz_min_2 - 1, nz_part_2 + kz
	DO k = ny_min_2 - 1, ny_part_2 + ky
		DO j = nx_min_2 - 1, nx_part_2 + kx

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
!$ACC END PARALLEL
!$OMP END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get flux !
IF(dual_energy) THEN 
	!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
	!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT)
	DO l = nz_min_2 - 1, nz_part_2 + kz
		DO k = ny_min_2 - 1, ny_part_2 + ky
			DO j = nx_min_2 - 1, nx_part_2 + kx
				uL2 (ieps2,j,k,l) = primL2(ieps2,j,k,l)
				uR2 (ieps2,j,k,l) = primR2(ieps2,j,k,l)
				fluxL2 (ieps2,j,k,l) = uL2 (ieps2,j,k,l) * primL2(ivn,j,k,l)
				fluxR2 (ieps2,j,k,l) = uR2 (ieps2,j,k,l) * primR2(ivn,j,k,l)
				fluxL2 (ieps2,j,k,l) = fluxL2 (ieps2,j,k,l) + primL2(itau2,j,k,l) * primL2(ivn,j,k,l)
				fluxR2 (ieps2,j,k,l) = fluxR2 (ieps2,j,k,l) + primR2(itau2,j,k,l) * primR2(ivn,j,k,l)
    	END DO
  	END DO
	END DO
	!$ACC END PARALLEL
	!$OMP END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL system_clock(time_end)
#ifdef DEBUG
WRITE(*,*) 'buildstates = ', REAL(time_end - time_start) / rate
#endif

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

! REAL !
REAL*8 :: dflux_2

! Integer !
INTEGER :: i, j, k, l

! Check timing with or without openmp
INTEGER :: time_start, time_end
INTEGER :: cr
REAL*8 :: rate
CALL system_clock(count_rate=cr)
rate = REAL(cr)

CALL system_clock(time_start)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now calculate flux gradient, choose according to the coordinate system !

!$OMP PARALLEL PRIVATE(dflux_2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! x-direciton !
IF(dir_in == x_dir) THEN
	IF(coordinate_flag == 0) THEN
		!$OMP DO COLLAPSE(4) SCHEDULE(STATIC) 
		!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT) PRIVATE(dflux_2)
		DO l = nz_min_2, nz_part_2
			DO k = ny_min_2, ny_part_2
				DO j = nx_min_2, nx_part_2
					DO i = imin2, ibx - 1
						dflux_2 = (flux_2 (i,j,k,l) - flux_2 (i,j-1,k,l)) / (dx2(j))
						l2(i,j,k,l) = - dflux_2
					END DO
				END DO
			END DO
		END DO
		!$ACC END PARALLEL
		!$OMP END DO
	ELSEIF(coordinate_flag == 1) THEN
		!$OMP DO COLLAPSE(4) SCHEDULE(STATIC)
		!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT) PRIVATE(dflux_2)
		DO l = nz_min_2, nz_part_2
			DO k = ny_min_2, ny_part_2
				DO j = nx_min_2, nx_part_2
					DO i = imin2, ibx - 1
						dflux_2 = (xF2(j)*flux_2 (i,j,k,l) - xF2(j-1)*flux_2 (i,j-1,k,l)) &
															/ (0.5D0*dx2_sq(j))
						l2(i,j,k,l) = - dflux_2
					END DO
				END DO
			END DO
		END DO
		!$ACC END PARALLEL
		!$OMP END DO
	ELSEIF(coordinate_flag == 2) THEN
		!$OMP DO COLLAPSE(4) SCHEDULE(STATIC)
		!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT) PRIVATE(dflux_2)
		DO l = nz_min_2, nz_part_2
			DO k = ny_min_2, ny_part_2
				DO j = nx_min_2, nx_part_2
					DO i = imin2, ibx - 1
						dflux_2 = (xF2(j)**2*flux_2 (i,j,k,l) - xF2(j-1)**2*flux_2 (i,j-1,k,l)) & 
															/ (dx2_cb(j)/3.0D0)
						l2(i,j,k,l) = - dflux_2
					END DO
				END DO
			END DO
		END DO
		!$ACC END PARALLEL
		!$OMP END DO
	END IF

! y-direciton !
ELSEIF(dir_in == y_dir) THEN
	IF(coordinate_flag == 0) THEN
		!$OMP DO COLLAPSE(4) SCHEDULE(STATIC)
		!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT) PRIVATE(dflux_2)
		DO l = nz_min_2, nz_part_2
			DO k = ny_min_2, ny_part_2
				DO j = nx_min_2, nx_part_2
					DO i = imin2, ibx - 1
						dflux_2 = (flux_2 (i,j,k,l) - flux_2 (i,j,k-1,l)) / (dy2(k))
						l2(i,j,k,l) = l2(i,j,k,l) - dflux_2
					END DO
				END DO
			END DO
		END DO
		!$ACC END PARALLEL
		!$OMP END DO
	ELSEIF(coordinate_flag == 1) THEN
		!$OMP DO COLLAPSE(4) SCHEDULE(STATIC)
		!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT) PRIVATE(dflux_2)
		DO l = nz_min_2, nz_part_2
			DO k = ny_min_2, ny_part_2
				DO j = nx_min_2, nx_part_2
					DO i = imin2, ibx - 1
						dflux_2 = (flux_2 (i,j,k,l) - flux_2 (i,j,k-1,l)) / (x2(j)*dy2(k))
						l2(i,j,k,l) = l2(i,j,k,l) - dflux_2
					END DO
				END DO
			END DO
		END DO
		!$ACC END PARALLEL
		!$OMP END DO
	ELSEIF(coordinate_flag == 2) THEN
		!$OMP DO COLLAPSE(4) SCHEDULE(STATIC)
		!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT) PRIVATE(dflux_2)
		DO l = nz_min_2, nz_part_2
			DO k = ny_min_2, ny_part_2
				DO j = nx_min_2, nx_part_2
					DO i = imin2, ibx - 1
						dflux_2 = (SIN(yF2(k))*flux_2 (i,j,k,l) - SIN(yF2(k-1))*flux_2 (i,j,k-1,l)) & 
															/ (x2bar(j)*dcos2(k))
						l2(i,j,k,l) = l2(i,j,k,l) - dflux_2
					END DO
				END DO
			END DO
		END DO
		!$ACC END PARALLEL
		!$OMP END DO
	END IF

! z-direciton !
ELSEIF(dir_in == z_dir) THEN
	IF(coordinate_flag == 0) THEN
		!$OMP DO COLLAPSE(4) SCHEDULE(STATIC)
		!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT) PRIVATE(dflux_2)
		DO l = nz_min_2, nz_part_2
			DO k = ny_min_2, ny_part_2
				DO j = nx_min_2, nx_part_2
					DO i = imin2, ibx - 1
						dflux_2 = (flux_2 (i,j,k,l) - flux_2 (i,j,k,l-1)) / dz2(l)
						l2(i,j,k,l) = l2(i,j,k,l) - dflux_2
					END DO
				END DO
			END DO
		END DO
		!$ACC END PARALLEL
		!$OMP END DO
	ELSEIF(coordinate_flag == 1) THEN
		!$OMP DO COLLAPSE(4) SCHEDULE(STATIC)
		!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT) PRIVATE(dflux_2)
		DO l = nz_min_2, nz_part_2
			DO k = ny_min_2, ny_part_2
				DO j = nx_min_2, nx_part_2
					DO i = imin2, ibx - 1
						dflux_2 = (flux_2 (i,j,k,l) - flux_2 (i,j,k,l-1)) / dz2(l)
						l2(i,j,k,l) = l2(i,j,k,l) - dflux_2
					END DO
				END DO
			END DO
		END DO
		!$ACC END PARALLEL
		!$OMP END DO
	ELSEIF(coordinate_flag == 2) THEN
		!$OMP DO COLLAPSE(4) SCHEDULE(STATIC)
		!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT) PRIVATE(dflux_2)
		DO l = nz_min_2, nz_part_2
			DO k = ny_min_2, ny_part_2
				DO j = nx_min_2, nx_part_2
					DO i = imin2, ibx - 1
						dflux_2 = (flux_2 (i,j,k,l) - flux_2 (i,j,k,l-1)) & 
															* dy2(k) / (x2bar(j)*dz2(l)*dcos2(k)) 
						l2(i,j,k,l) = l2(i,j,k,l) - dflux_2
					END DO
				END DO
			END DO
		END DO
		!$ACC END PARALLEL
		!$OMP END DO
	END IF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP END PARALLEL 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL system_clock(time_end)
#ifdef DEBUG
WRITE(*,*) 'flux_diff = ', REAL(time_end - time_start) / rate
#endif

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
REAL*8 :: dbxdt_l, dbxdt_r
REAL*8 :: dbydt_l, dbydt_r
REAL*8 :: dbzdt_l, dbzdt_r

! Dummy !
REAL*8 :: dbdx, dbdy, dbdz, divb

! Check timing with or without openmp
INTEGER :: time_start, time_end
INTEGER :: cr
REAL*8 :: rate

CALL system_clock(count_rate=cr)
rate = REAL(cr)

CALL system_clock(time_start)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP PARALLEL PRIVATE(dbxdt_l, dbxdt_r, dbydt_l, dbydt_r, dbzdt_l, dbzdt_r, dbdx, dbdy, dbdz, divb) REDUCTION(MAX:maxDivB)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Final step, get rungekutta operator, LHS of the hydro equation !
!$OMP DO COLLAPSE(4) SCHEDULE(STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT)
DO l = nz_min_2, nz_part_2
	DO k = ny_min_2, ny_part_2
		DO j = nx_min_2, nx_part_2
			DO i = imin2, ibx - 1
				l2(i,j,k,l) = l2(i,j,k,l) + sc2(i,j,k,l)
			END DO
		END DO
	END DO
END DO
!$ACC END PARALLEL
!$OMP END DO

! Initialize !
maxDivB = 0.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Magnetic field, special treatment for cylindrical coordinates !
IF(coordinate_flag == 0) THEN
	!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
	!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) &
	!$ACC PRIVATE(dbxdt_l, dbxdt_r, dbydt_l, dbydt_r, dbzdt_l, dbzdt_r, dbdx, dbdy, dbdz, divb) REDUCTION(MAX:maxDivB)
	DO l = nz_min_2, nz_part_2
		DO k = ny_min_2, ny_part_2
			DO j = nx_min_2, nx_part_2

				! dbx/dt !
				dbxdt_r = - ((efield_z(j,k,l) - efield_z(j,k-1,l))/(dy2(k)) - (efield_y(j,k,l) - efield_y(j,k,l-1))/(dz2(l)))
				dbxdt_l = - ((efield_z(j-1,k,l) - efield_z(j-1,k-1,l))/(dy2(k)) - (efield_y(j-1,k,l) - efield_y(j-1,k,l-1))/(dz2(l)))
				l2(ibx,j,k,l) = 0.5D0*(dbxdt_l + dbxdt_r)

				! dby/dt !
				dbydt_r = - ((efield_x(j,k,l) - efield_x(j,k,l-1))/(dz2(l)) - (efield_z(j,k,l) - efield_z(j-1,k,l))/(dx2(j)))
				dbydt_l = - ((efield_x(j,k-1,l) - efield_x(j,k-1,l-1))/(dz2(l)) - (efield_z(j,k-1,l) - efield_z(j-1,k-1,l))/(dx2(j)))
				l2(iby,j,k,l) = 0.5D0*(dbydt_l + dbydt_r)

				! dbz/dt !
				dbzdt_r = - ((efield_y(j,k,l) - efield_y(j-1,k,l))/(dx2(j)) - (efield_x(j,k,l) - efield_x(j,k-1,l))/(dy2(k)))
				dbzdt_l = - ((efield_y(j,k,l-1) - efield_y(j-1,k,l-1))/(dx2(j)) - (efield_x(j,k,l-1) - efield_x(j,k-1,l-1))/(dy2(k)))
				l2(ibz,j,k,l) = 0.5D0*(dbzdt_l + dbzdt_r)

				! divergence B !
        dbdx = (dbxdt_r - dbxdt_l)/(dx2(j))
        dbdy = (dbydt_r - dbydt_l)/(dy2(k))
        dbdz = (dbzdt_r - dbzdt_l)/(dz2(l))
        divb = dbdx + dbdy + dbdz
        maxDivB = max(maxDivB, divb)

			END DO
		END DO
	END DO
	!$ACC END PARALLEL
	!$OMP END DO
ELSEIF(coordinate_flag == 1) THEN
	!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
	!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) &
	!$ACC PRIVATE(dbxdt_l, dbxdt_r, dbydt_l, dbydt_r, dbzdt_l, dbzdt_r, dbdx, dbdy, dbdz, divb) REDUCTION(MAX:maxDivB)
	DO l = nz_min_2, nz_part_2
		DO k = ny_min_2, ny_part_2
			DO j = nx_min_2, nx_part_2

				! dbx/dt !
				dbxdt_r = - ((efield_z(j,k,l) - efield_z(j,k-1,l))/(xF2(j)*dy2(k)) - (efield_y(j,k,l) - efield_y(j,k,l-1))/(dz2(l)))
				dbxdt_l = - ((efield_z(j-1,k,l) - efield_z(j-1,k-1,l))/(xF2(j-1)*dy2(k)) - (efield_y(j-1,k,l) - efield_y(j-1,k,l-1))/(dz2(l)))
				l2(ibx,j,k,l) = (xF2(j-1)*dbxdt_l + xF2(j)*dbxdt_r)/(xF2(j) + xF2(j-1))

				! dby/dt !
				dbydt_r = - ((efield_x(j,k,l) - efield_x(j,k,l-1))/(dz2(l)) - (efield_z(j,k,l) - efield_z(j-1,k,l))/(dx2(j)))
				dbydt_l = - ((efield_x(j,k-1,l) - efield_x(j,k-1,l-1))/(dz2(l)) - (efield_z(j,k-1,l) - efield_z(j-1,k-1,l))/(dx2(j)))
				l2(iby,j,k,l) = 0.5D0*(dbydt_l + dbydt_r)

				! dbz/dt !
				dbzdt_r = - ((xF2(j)*efield_y(j,k,l) - xF2(j-1)*efield_y(j-1,k,l))/(dx2(j)) - (efield_x(j,k,l) - efield_x(j,k-1,l))/(dy2(k)))/x2(j)
				dbzdt_l = - ((xF2(j)*efield_y(j,k,l-1) - xF2(j-1)*efield_y(j-1,k,l-1))/(dx2(j)) - (efield_x(j,k,l-1) - efield_x(j,k-1,l-1))/(dy2(k)))/x2(j)
				l2(ibz,j,k,l) = 0.5D0*(dbzdt_l + dbzdt_r)

				! divegence b !
        dbdx = (xF2(j)*dbxdt_r - xF2(j-1)*dbxdt_l)/(x2(j)*dx2(j))
        dbdy = (dbydt_r - dbydt_l)/(x2(j)*dy2(k))
        dbdz = (dbzdt_r - dbzdt_l)/(dz2(l))
        divb = dbdx + dbdy + dbdz
        maxDivB = max(maxDivB, divb)

			END DO
		END DO
	END DO
	!$ACC END PARALLEL
	!$OMP END DO
ELSEIF(coordinate_flag == 2) THEN
	!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
	!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) &
	!$ACC PRIVATE(dbxdt_l, dbxdt_r, dbydt_l, dbydt_r, dbzdt_l, dbzdt_r, dbdx, dbdy, dbdz, divb) REDUCTION(MAX:maxDivB)
	DO l = nz_min_2, nz_part_2
		DO k = ny_min_2, ny_part_2
			DO j = nx_min_2, nx_part_2

				! dbx/dt !
				dbxdt_r = - ((SIN(yF2(k))*efield_z(j,k,l) - SIN(yF2(k-1))*efield_z(j,k-1,l))/(dcos2(k)) - (efield_y(j,k,l) - efield_y(j,k,l-1))/(dz2(l))*SIN(y2(k)))/(xF2(j))
				dbxdt_l = - ((SIN(yF2(k))*efield_z(j-1,k,l) - SIN(yF2(k-1))*efield_z(j-1,k-1,l))/(dcos2(k)) - (efield_y(j-1,k,l) - efield_y(j-1,k,l-1))/(dz2(l))*SIN(y2(k)))/(xF2(j-1))
				l2(ibx,j,k,l) = 1.5D0*dx2(j)/dx2_cb(j)*(xF2(j-1)**2*dbxdt_l + xF2(j)**2*dbxdt_r)

				! dby/dt !
				dbydt_r = - ((efield_x(j,k,l) - efield_x(j,k,l-1))/(SIN(yF2(k))*dz2(l)) - (xF2(j)*efield_z(j,k,l) - xF2(j-1)*efield_z(j-1,k,l))/(dx2(j)))/x2(j)
				dbydt_l = - ((efield_x(j,k-1,l) - efield_x(j,k-1,l-1))/(SIN(yF2(k-1))*dz2(l)) - (xF2(j)*efield_z(j,k-1,l) - xF2(j-1)*efield_z(j-1,k-1,l))/(dx2(j)))/x2(j)
				l2(iby,j,k,l) = 0.5D0*dy2(k)/dcos2(k)*(SIN(yF2(k-1))*dbydt_l + SIN(yF2(k))*dbydt_r)

				! dbz/dt !
				dbzdt_r = - ((xF2(j)*efield_y(j,k,l) - xF2(j-1)*efield_y(j-1,k,l))/(dx2(j)) - (efield_x(j,k,l) - efield_x(j,k-1,l))/(dy2(k)))/x2(j)
				dbzdt_l = - ((xF2(j)*efield_y(j,k,l-1) - xF2(j-1)*efield_y(j-1,k,l-1))/(dx2(j)) - (efield_x(j,k,l-1) - efield_x(j,k-1,l-1))/(dy2(k)))/x2(j)
				l2(ibz,j,k,l) = 0.5D0*(dbzdt_l + dbzdt_r)

				! divergence b !
        dbdx = (xF2(j)**2*dbxdt_r - xF2(j-1)**2*dbxdt_l)/(x2(j)**2*dx2(j))
        dbdy = (SIN(yF2(k))*dbydt_r - SIN(yF2(k-1))*dbydt_l)/(x2(j)*dcos2(k))!*SIN(y2(k))*dy2(k))
        dbdz = (dbzdt_r - dbzdt_l)/(x2(j)*SIN(y2(k))*dz2(l))
        divb = dbdx + dbdy + dbdz
        maxDivB = max(maxDivB, divb)

			END DO
		END DO
	END DO
	!$ACC END PARALLEL
	!$OMP END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL system_clock(time_end)
#ifdef DEBUG
WRITE(*,*) 'l_operator = ', REAL(time_end - time_start) / rate
#endif

END SUBROUTINE