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
USE TVD_MODULE
USE PPMC_MODULE
USE WENO_MODULE
USE RIEMANN_MODULE
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find source term !

! Predefined source term !
CALL GET_SOURCE

! Custom source term !
CALL CUSTOM_SOURCE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do the x-sweep !

! Reconstruction !
IF(tvdmm_flag) THEN
	CALL tvd_reconx(minmod)
ELSEIF(tvdvl_flag) THEN
	CALL tvd_reconx(vanleer)
ELSEIF(tvdmc_flag) THEN
	CALL tvd_reconx(mcentral)
ELSEIF(ppmc_flag) THEN
	CALL ppmc_reconx
ELSEIF(weno_flag) THEN
	CALL weno_reconx
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

! Back up fluxes !
CALL GETFLUX_X

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do the y-sweep !
IF(n_dim > 1) THEN

	! Reconstruction !
	IF(tvdmm_flag) THEN
		CALL tvd_recony(minmod)
	ELSEIF(tvdvl_flag) THEN
		CALL tvd_recony(vanleer)
	ELSEIF(tvdmc_flag) THEN
		CALL tvd_recony(mcentral)
	ELSEIF(ppmc_flag) THEN
		CALL ppmc_recony
	ELSEIF(weno_flag) THEN
		CALL weno_recony
	END IF

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

	! Back up fluxes !
	CALL GETFLUX_Y

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do the z-sweep !

IF(n_dim > 2) THEN

	! Reconstruction !
	IF(tvdmm_flag) THEN
		CALL tvd_reconz(minmod)
	ELSEIF(tvdvl_flag) THEN
		CALL tvd_reconz(vanleer)
	ELSEIF(tvdmc_flag) THEN
		CALL tvd_reconz(mcentral)
	ELSEIF(ppmc_flag) THEN
		CALL ppmc_reconz
	ELSEIF(weno_flag) THEN
		CALL weno_reconz
	END IF

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

	! Back up fluxes !
	CALL GETFLUX_Z

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
USE MHD_MODULE
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k, l

! Pressure gradients !
REAL*8 :: bsquare

! Check timing with or without openmp
#ifdef DEBUG
INTEGER :: time_start, time_end
INTEGER :: cr
REAL*8 :: rate
CALL system_clock(count_rate=cr)
rate = REAL(cr)
CALL system_clock(time_start)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now do the normal matter

! Geometric sources
!$OMP PARALLEL PRIVATE(bsquare)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize !
!$OMP DO COLLAPSE(4) SCHEDULE(STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT) 
DO l = 1, nz
	DO k = 1, ny
		DO j = 1, nx
			DO i = imin, imax
				sc(i,j,k,l) = 0.0D0
			END DO
		END DO
	END DO
END DO
!$ACC END PARALLEL
!$OMP END DO

! Choose coordinate system, add geometric source term !
IF(coordinate_flag == 1) THEN
	!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
	!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) PRIVATE(bsquare)
	DO l = 1, nz
		DO k = 1, ny
			DO j = 1, nx
				bsquare = dot_product(bcell(ibx:ibz,j,k,l),bcell(ibx:ibz,j,k,l))
				sc(ivx,j,k,l) = sc(ivx,j,k,l) + (prim(itau,j,k,l) + prim(irho,j,k,l)*prim(ivy,j,k,l)*prim(ivy,j,k,l) & 
																									+ 0.5D0*bsquare - bcell(iby,j,k,l)*bcell(iby,j,k,l))/(x(j))
																																											
				sc(ivy,j,k,l) = sc(ivy,j,k,l) - (prim(irho,j,k,l)*prim(ivx,j,k,l)*prim(ivy,j,k,l) & 
																									- bcell(ibx,j,k,l)*bcell(iby,j,k,l))/(x(j))
			END DO
		END DO
	END DO
	!$ACC END PARALLEL
	!$OMP END DO
ELSEIF(coordinate_flag == 2) THEN
	!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
	!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) PRIVATE(bsquare)
	DO l = 1, nz
		DO k = 1, ny
			DO j = 1, nx
				bsquare = dot_product(bcell(ibx:ibz,j,k,l),bcell(ibx:ibz,j,k,l))
				sc(ivx,j,k,l) = sc(ivx,j,k,l) + (2.0D0*prim(itau,j,k,l) + bcell(ibx,j,k,l)*bcell(ibx,j,k,l) &
																								+ prim(irho,j,k,l)*(prim(ivy,j,k,l)*prim(ivy,j,k,l) &
																								+ prim(ivz,j,k,l)*prim(ivz,j,k,l))) & 
																								* (3.0D0*x(j)*dx(j)/dx_cb(j)) 

				sc(ivy,j,k,l) = sc(ivy,j,k,l) + (prim(itau,j,k,l) + prim(irho,j,k,l)*prim(ivz,j,k,l) &
																								* prim(ivz,j,k,l) + 0.5D0*bsquare - bcell(ibz,j,k,l)*bcell(ibz,j,k,l)) &
																								* (3.0D0*x(j)*dx(j)/dx_cb(j))*(dsine(k)/dcose(k)) & 
																						
																								- (prim(irho,j,k,l)*prim(ivx,j,k,l)*prim(ivy,j,k,l) & 
																								- bcell(ibx,j,k,l)*bcell(iby,j,k,l))*(3.0D0*x(j)*dx(j)/dx_cb(j)) 

				sc(ivz,j,k,l) = sc(ivz,j,k,l) - (prim(irho,j,k,l)*prim(ivx,j,k,l)*prim(ivz,j,k,l) & 
																								- bcell(ibx,j,k,l)*bcell(ibz,j,k,l))*(3.0D0*x(j)*dx(j)/dx_cb(j)) & 

																								- (prim(irho,j,k,l)*prim(ivy,j,k,l)*prim(ivz,j,k,l) & 
																								- bcell(iby,j,k,l)*bcell(ibz,j,k,l))*(3.0D0*x(j)*dx(j)/dx_cb(j))*(dsine(k)/dcose(k))
			END DO
		END DO
	END DO
	!$ACC END PARALLEL
	!$OMP END DO
END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$OMP END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'getsource = ', REAL(time_end - time_start) / rate
#endif

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine built the states for left and right edges for the horizontal directions !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BUILDSTATES(dir_in)
USE DEFINITION 
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
#ifdef DEBUG
INTEGER :: time_start, time_end
INTEGER :: cr
REAL*8 :: rate
CALL system_clock(count_rate=cr)
rate = REAL(cr)
CALL system_clock(time_start)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign !
kx = 1
ky = 1
kz = 1

! Assign !
IF(dir_in == x_dir) THEN
	ivn = ivx
	ivt1 = ivy
	ivt2 = ivz
	ibn = ibx
	ibt1 = iby
	ibt2 = ibz
	kx = 0
ELSEIF(dir_in == y_dir) THEN
	ivn = ivy
	ivt1 = ivz
	ivt2 = ivx
	ibn = iby
	ibt1 = ibz
	ibt2 = ibx
	ky = 0
ELSEIF(dir_in == z_dir) THEN
	ivn = ivz
	ivt1 = ivx
	ivt2 = ivy
	ibn = ibz
	ibt1 = ibx
	ibt2 = iby
	kz = 0
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get flux !

!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(v2L, v2R, b2L, b2R, vbL, vbR)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE (3) DEFAULT(PRESENT) PRIVATE(v2L, v2R, b2L, b2R, vbL, vbR)
DO l = 0, nz + kz
	DO k = 0, ny + ky
		DO j = 0, nx + kx

				! get dot product !
				b2L = dot_product(primL(ibx:ibz,j,k,l), primL(ibx:ibz,j,k,l))
				b2R = dot_product(primR(ibx:ibz,j,k,l), primR(ibx:ibz,j,k,l))
				v2L = dot_product(primL(ivx:ivz,j,k,l), primL(ivx:ivz,j,k,l))
				v2R = dot_product(primR(ivx:ivz,j,k,l), primR(ivx:ivz,j,k,l))
				vbL = dot_product(primL(ibx:ibz,j,k,l), primL(ivx:ivz,j,k,l))
				vbR = dot_product(primR(ibx:ibz,j,k,l), primR(ivx:ivz,j,k,l))

				! conservative variables !
				uL (imin:ibx-1,j,k,l) = primL(imin:ibx-1,j,k,l)*primL(irho,j,k,l)
				uR (imin:ibx-1,j,k,l) = primR(imin:ibx-1,j,k,l)*primR(irho,j,k,l)
				uL (irho,j,k,l) = primL(irho,j,k,l)
				uR (irho,j,k,l) = primR(irho,j,k,l)

				! Magnetic fields !
				uL (ibx:ibz,j,k,l) = primL(ibx:ibz,j,k,l)
				uR (ibx:ibz,j,k,l) = primR(ibx:ibz,j,k,l)

				! NM energy equation !
				uL (itau,j,k,l) = primL(irho,j,k,l)*(0.5D0*v2L + epsL(j,k,l)) + 0.5D0*b2L
				uR (itau,j,k,l) = primR(irho,j,k,l)*(0.5D0*v2R + epsR(j,k,l)) + 0.5D0*b2R

				! For NM flux !
				fluxL (imin:ibx-1,j,k,l) = uL (imin:ibx-1,j,k,l) * primL(ivn,j,k,l)
				fluxR (imin:ibx-1,j,k,l) = uR (imin:ibx-1,j,k,l) * primR(ivn,j,k,l)

				! Add the pressure term to x-momentum equation !
				fluxL (ivn,j,k,l) = fluxL (ivn,j,k,l) + primL(itau,j,k,l) + 0.5D0*(b2L) - primL(ibn,j,k,l)*primL(ibn,j,k,l)
				fluxR (ivn,j,k,l) = fluxR (ivn,j,k,l) + primR(itau,j,k,l) + 0.5D0*(b2R) - primR(ibn,j,k,l)*primR(ibn,j,k,l)

				! Other momentum fluxes !
				fluxL(ivt1,j,k,l) = fluxL(ivt1,j,k,l) - primL(ibn,j,k,l)*primL(ibt1,j,k,l)
				fluxR(ivt1,j,k,l) = fluxR(ivt1,j,k,l) - primR(ibn,j,k,l)*primR(ibt1,j,k,l)
				fluxL(ivt2,j,k,l) = fluxL(ivt2,j,k,l) - primL(ibn,j,k,l)*primL(ibt2,j,k,l)
				fluxR(ivt2,j,k,l) = fluxR(ivt2,j,k,l) - primR(ibn,j,k,l)*primR(ibt2,j,k,l)

				! Add the presusre work done term to the energy equation             
				fluxL (itau,j,k,l) = fluxL (itau,j,k,l) + (primL(itau,j,k,l) + 0.5D0*(b2L)) * primL(ivn,j,k,l) - primL(ibn,j,k,l)*vbL
				fluxR (itau,j,k,l) = fluxR (itau,j,k,l) + (primR(itau,j,k,l) + 0.5D0*(b2R)) * primR(ivn,j,k,l) - primR(ibn,j,k,l)*vbR

				! magentic field fluxes !
				fluxL(ibn,j,k,l) = 0.0D0
				fluxR(ibn,j,k,l) = 0.0D0
				fluxL(ibt1,j,k,l) = (primL(ivn,j,k,l)*primL(ibt1,j,k,l) - primL(ivt1,j,k,l)*primL(ibn,j,k,l))
				fluxR(ibt1,j,k,l) = (primR(ivn,j,k,l)*primR(ibt1,j,k,l) - primR(ivt1,j,k,l)*primR(ibn,j,k,l))
				fluxL(ibt2,j,k,l) = (primL(ivn,j,k,l)*primL(ibt2,j,k,l) - primL(ivt2,j,k,l)*primL(ibn,j,k,l))
				fluxR(ibt2,j,k,l) = (primR(ivn,j,k,l)*primR(ibt2,j,k,l) - primR(ivt2,j,k,l)*primR(ibn,j,k,l))

		END DO
	END DO
END DO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'buildstates = ', REAL(time_end - time_start) / rate
#endif

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine built the states for left and right edges for the horizontal directions !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FLUX_DIFF(dir_in)
USE DEFINITION 
IMPLICIT NONE

! Integer !
INTEGER, INTENT(IN) :: dir_in

! REAL !
REAL*8 :: dflux

! Integer !
INTEGER :: i, j, k, l

! Check timing with or without openmp
#ifdef DEBUG
INTEGER :: time_start, time_end
INTEGER :: cr
REAL*8 :: rate
CALL system_clock(count_rate=cr)
rate = REAL(cr)
CALL system_clock(time_start)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now calculate flux gradient, choose according to the coordinate system !

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! x-direciton !
IF(dir_in == x_dir) THEN
	IF(coordinate_flag == 0) THEN
		!$OMP PARALLEL DO COLLAPSE(4) SCHEDULE(STATIC) PRIVATE(dflux)
		!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT) PRIVATE(dflux)
		DO l = 1, nz
			DO k = 1, ny
				DO j = 1, nx
					DO i = imin, ibx - 1
						dflux = (flux (i,j,k,l) - flux (i,j-1,k,l)) / (dx(j))
						l_rk(i,j,k,l) = - dflux
					END DO
				END DO
			END DO
		END DO
		!$ACC END PARALLEL
		!$OMP END PARALLEL DO
	ELSEIF(coordinate_flag == 1) THEN
		!$OMP PARALLEL DO COLLAPSE(4) SCHEDULE(STATIC) PRIVATE(dflux)
		!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT) PRIVATE(dflux)
		DO l = 1, nz
			DO k = 1, ny
				DO j = 1, nx
					DO i = imin, ibx - 1
						dflux = (xF(j)*flux (i,j,k,l) - xF(j-1)*flux (i,j-1,k,l)) / (x(j)*dx(j))
						l_rk(i,j,k,l) = - dflux
					END DO
				END DO
			END DO
		END DO
		!$ACC END PARALLEL
		!$OMP END PARALLEL DO
	ELSEIF(coordinate_flag == 2) THEN
		!$OMP PARALLEL DO COLLAPSE(4) SCHEDULE(STATIC) PRIVATE(dflux)
		!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT) PRIVATE(dflux)
		DO l = 1, nz
			DO k = 1, ny
				DO j = 1, nx
					DO i = imin, ibx - 1
						dflux = (xF(j)*xF(j)*flux (i,j,k,l) - xF(j-1)*xF(j-1)*flux (i,j-1,k,l)) / (dx_cb(j)/3.0D0)
						l_rk(i,j,k,l) = - dflux
					END DO
				END DO
			END DO
		END DO
		!$ACC END PARALLEL
		!$OMP END PARALLEL DO
	END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! y-direciton !
ELSEIF(dir_in == y_dir) THEN
	IF(coordinate_flag == 0) THEN
		!$OMP PARALLEL DO COLLAPSE(4) SCHEDULE(STATIC) PRIVATE(dflux)
		!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT) PRIVATE(dflux)
		DO l = 1, nz
			DO k = 1, ny
				DO j = 1, nx
					DO i = imin, ibx - 1
						dflux = (flux (i,j,k,l) - flux (i,j,k-1,l)) / (dy(k))
						l_rk(i,j,k,l) = l_rk(i,j,k,l) - dflux
					END DO
				END DO
			END DO
		END DO
		!$ACC END PARALLEL
		!$OMP END PARALLEL DO
	ELSEIF(coordinate_flag == 1) THEN
		!$OMP PARALLEL DO COLLAPSE(4) SCHEDULE(STATIC) PRIVATE(dflux)
		!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT) PRIVATE(dflux)
		DO l = 1, nz
			DO k = 1, ny
				DO j = 1, nx
					DO i = imin, ibx - 1
						dflux = (flux (i,j,k,l) - flux (i,j,k-1,l)) / (x(j)*dy(k))
						l_rk(i,j,k,l) = l_rk(i,j,k,l) - dflux
					END DO
				END DO
			END DO
		END DO
		!$ACC END PARALLEL
		!$OMP END PARALLEL DO
	ELSEIF(coordinate_flag == 2) THEN
		!$OMP PARALLEL DO COLLAPSE(4) SCHEDULE(STATIC) PRIVATE(dflux)
		!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT) PRIVATE(dflux)
		DO l = 1, nz
			DO k = 1, ny
				DO j = 1, nx
					DO i = imin, ibx - 1
						dflux = (sinf(k)*flux (i,j,k,l) - sinf(k-1)*flux (i,j,k-1,l)) / (xbar(j)*dcose(k))
						l_rk(i,j,k,l) = l_rk(i,j,k,l) - dflux
					END DO
				END DO
			END DO
		END DO
		!$ACC END PARALLEL
		!$OMP END PARALLEL DO
	END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! z-direciton !
ELSEIF(dir_in == z_dir) THEN
	IF(coordinate_flag == 0) THEN
		!$OMP PARALLEL DO COLLAPSE(4) SCHEDULE(STATIC) PRIVATE(dflux)
		!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT) PRIVATE(dflux)
		DO l = 1, nz
			DO k = 1, ny
				DO j = 1, nx
					DO i = imin, ibx - 1
						dflux = (flux (i,j,k,l) - flux (i,j,k,l-1)) / dz(l)
						l_rk(i,j,k,l) = l_rk(i,j,k,l) - dflux
					END DO
				END DO
			END DO
		END DO
		!$ACC END PARALLEL
		!$OMP END PARALLEL DO
	ELSEIF(coordinate_flag == 1) THEN
		!$OMP PARALLEL DO COLLAPSE(4) SCHEDULE(STATIC) PRIVATE(dflux)
		!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT) PRIVATE(dflux)
		DO l = 1, nz
			DO k = 1, ny
				DO j = 1, nx
					DO i = imin, ibx - 1
						dflux = (flux (i,j,k,l) - flux (i,j,k,l-1)) / dz(l)
						l_rk(i,j,k,l) = l_rk(i,j,k,l) - dflux
					END DO
				END DO
			END DO
		END DO
		!$ACC END PARALLEL
		!$OMP END PARALLEL DO
	ELSEIF(coordinate_flag == 2) THEN
		!$OMP PARALLEL DO COLLAPSE(4) SCHEDULE(STATIC) PRIVATE(dflux)
		!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT) PRIVATE(dflux)
		DO l = 1, nz
			DO k = 1, ny
				DO j = 1, nx
					DO i = imin, ibx - 1
						dflux = (flux (i,j,k,l) - flux (i,j,k,l-1)) * dy(k) / (xbar(j)*dz(l)*dcose(k)) 
						l_rk(i,j,k,l) = l_rk(i,j,k,l) - dflux
					END DO
				END DO
			END DO
		END DO
		!$ACC END PARALLEL
		!$OMP END PARALLEL DO
	END IF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'flux_diff = ', REAL(time_end - time_start) / rate
#endif

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine built the states for left and right edges for the horizontal directions !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE l_operator
USE DEFINITION 
USE MHD_MODULE
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k, l, lc, lcm1

! Real !
REAL*8 :: rbar

! Check timing with or without openmp
#ifdef DEBUG
INTEGER :: time_start, time_end
INTEGER :: cr
REAL*8 :: rate
CALL system_clock(count_rate=cr)
rate = REAL(cr)
CALL system_clock(time_start)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP PARALLEL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Final step, get rungekutta operator, LHS of the hydro equation !
!$OMP DO COLLAPSE(4) SCHEDULE(STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT)
DO l = 1, nz
	DO k = 1, ny
		DO j = 1, nx
			DO i = imin, ibx - 1
				l_rk(i,j,k,l) = l_rk(i,j,k,l) + sc(i,j,k,l)
			END DO
		END DO
	END DO
END DO
!$ACC END PARALLEL
!$OMP END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Magnetic field, we do it by line-integral but not direct numerical derivative !
IF(coordinate_flag == 0) THEN
	!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
	!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) 
	DO l = 0, nz
		DO k = 0, ny
			DO j = 0, nx

				! dbx/dt !
				l_rk(ibx,j,k,l) = - ((efield_z(j,k,l) - efield_z(j,k-1,l))/(dy(k)) - (efield_y(j,k,l) - efield_y(j,k,l-1))/(dz(l)))

				! dby/dt !
				l_rk(iby,j,k,l) = - ((efield_x(j,k,l) - efield_x(j,k,l-1))/(dz(l)) - (efield_z(j,k,l) - efield_z(j-1,k,l))/(dx(j)))

				! dbz/dt !
				l_rk(ibz,j,k,l) = - ((efield_y(j,k,l) - efield_y(j-1,k,l))/(dx(j)) - (efield_x(j,k,l) - efield_x(j,k-1,l))/(dy(k)))

			END DO
		END DO
	END DO
	!$ACC END PARALLEL
	!$OMP END DO
ELSEIF(coordinate_flag == 1) THEN
	!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
	!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) 
	DO l = 0, nz
		DO k = 0, ny
			DO j = 0, nx

				! dbx/dt !
				l_rk(ibx,j,k,l) = - (efield_z(j,k,l) - efield_z(j,k-1,l))/(xF(j)*dy(k)) + (efield_y(j,k,l) - efield_y(j,k,l-1))/(dz(l))

				! dby/dt !
				l_rk(iby,j,k,l) = - (efield_x(j,k,l) - efield_x(j,k,l-1))/(dz(l)) + (efield_z(j,k,l) - efield_z(j-1,k,l))/(dx(j))

				! dbz/dt !
				l_rk(ibz,j,k,l) = - (xF(j)*efield_y(j,k,l) - xF(j-1)*efield_y(j-1,k,l))/(x(j)*dx(j)) + (efield_x(j,k,l) - efield_x(j,k-1,l))/(x(j)*dy(k))

			END DO
		END DO
	END DO
	!$ACC END PARALLEL
	!$OMP END DO
ELSEIF(coordinate_flag == 2) THEN
	!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
	!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT)
	DO l = 0, nz
		DO k = 0, ny
			DO j = 0, nx

				! dbx/dt !
				l_rk(ibx,j,k,l) = (- (sinf(k)*efield_z(j,k,l) - sinf(k-1)*efield_z(j,k-1,l)) + (efield_y(j,k,l) - efield_y(j,k,l-1))*(dy(k)/dz(l)))/(xF(j)*dcose(k)+small_num)
				
				! dby/dt !
				l_rk(iby,j,k,l) = (- (efield_x(j,k,l) - efield_x(j,k,l-1))/(xbar(j)*sinf(k)*dz(l)+small_num) + (xF(j)*efield_z(j,k,l) - xF(j-1)*efield_z(j-1,k,l))/(x(j)*dx(j)))

				! dbz/dt !
				l_rk(ibz,j,k,l) = (- (xF(j)*efield_y(j,k,l) - xF(j-1)*efield_y(j-1,k,l))/(x(j)*dx(j)) + (efield_x(j,k,l) - efield_x(j,k-1,l))/(xbar(j)*dy(k)))

			END DO
		END DO
	END DO
	!$ACC END PARALLEL
	!$OMP END DO
END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$OMP END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'l_operator = ', REAL(time_end - time_start) / rate
#endif

END SUBROUTINE
