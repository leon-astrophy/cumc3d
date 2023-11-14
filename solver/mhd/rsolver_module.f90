!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This files contain all the riemann solvers available for 
! simulating hydrodynamics				   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE RIEMANN_MODULE
USE DEFINITION
IMPLICIT NONE

! HLLC solver variables  
REAL*8, ALLOCATABLE, DIMENSION(:) :: u_hll, p_hll
REAL*8, ALLOCATABLE, DIMENSION(:) :: ustarL, ustarR
REAL*8, ALLOCATABLE, DIMENSION(:) :: usstarL, usstarR
REAL*8, ALLOCATABLE, DIMENSION(:) :: pstarL, pstarR
REAL*8, ALLOCATABLE, DIMENSION(:) :: psstarL, psstarR

! Left and right hydro-states for DM/NM !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: epsR, epsL

! Speed of sound !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: csL, csR

! Left and right fluxes, conserved quantity !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: primL, primR
REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: fluxL, fluxR, uL, uR

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! function for computing the signal speed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	REAL*8 function compute_signalspeed(cs, bn, bt1, bt2, rho)
	!$ACC ROUTINE SEQ
	implicit none
	REAL*8 :: cs, bn, bt1, bt2, rho
	REAL*8 :: a2_mhd, b2_mhd, a4_mhd, b4_mhd
	REAL*8 :: b2n_mhd, b2t1_mhd, b2t2_mhd
	a2_mhd = cs*cs
	a4_mhd = a2_mhd*a2_mhd
	b2n_mhd = (bn*bn/rho)
	b2t1_mhd = (bt1*bt1/rho)
	b2t2_mhd = (bt2*bt2/rho)
	b2_mhd = b2n_mhd + b2t1_mhd + b2t2_mhd
	b4_mhd = b2_mhd*b2_mhd
	compute_signalspeed = DSQRT(0.5D0*(a2_mhd + b2_mhd + & 
	                      DSQRT(MAX((a4_mhd + b4_mhd + 2.0d0*b2_mhd*a2_mhd) & 
												- 4.0D0*a2_mhd*b2n_mhd, 0.0d0))))
	end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! function for computing the HLL flux
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	REAL*8 function compute_hll(ul,ur,fl,fr,sl,sr)
	!$ACC ROUTINE SEQ
	implicit none
	REAL*8 :: ul, ur, fl, fr, sl, sr
	compute_hll = (sr*ur - sl*ul - (fr - fl))/(sr - sl)
	end function

	REAL*8 function compute_fluxhll(yl,yr,xl,xr,sl,sr)
	!$ACC ROUTINE SEQ
	implicit none
	REAL*8 :: yl, yr, xl, xr, sl, sr
	compute_fluxhll = (sr*yl-sl*yr+sl*sr*(xr-xl))/(sr-sl)
	end function

	REAL*8 function compute_roe(xl,xr,rhol,rhor)
	!$ACC ROUTINE SEQ
	implicit none
	REAL*8 :: xl,xr,rhol,rhor
	compute_roe = (DSQRT(rhol)*xl + DSQRT(rhor)*xr)/(DSQRT(rhol) + DSQRT(rhor))
	end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! function for computing the HLLC flux
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	REAL*8 function compute_plr(pl,pr,rhol,rhor,vl,vr,sl,sr,star)
	!$ACC ROUTINE SEQ
	implicit none
	REAL*8 :: pl,pr,rhol,rhor,vl,vr,sl,sr,star
	compute_plr = 0.5D0*(pl+pr+rhol*(sl-vl)*(star-vl)+rhor*(sr-vr)*(star-vr))
	end function

	REAL*8 function compute_sstar_hllc(pl,pr,rhol,rhor,vl,vr,sl,sr,bl,br)
	!$ACC ROUTINE SEQ
	implicit none
	REAL*8 :: pl, pr, rhol, rhor, vl, vr, sl, sr, bl, br
	compute_sstar_hllc = (rhor*vr*(sr - vr) - rhol*vl*(sl-vl) + & 
												pl - pr + br*br - bl*bl) / &
												(rhor*(sr-vr) - rhol*(sl-vl))
	end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Assign left and right states and fluxes for riemann problem !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BUILD_RIEMANN
USE DEFINITION 
IMPLICIT NONE

! HLLC variables !
ALLOCATE(ustarL(imin:imax))
ALLOCATE(ustarR(imin:imax))
ALLOCATE(pstarL(imin:imax))
ALLOCATE(pstarR(imin:imax))
ALLOCATE(u_hll(imin:imax))
ALLOCATE(p_hll(imin:imax))

! HLLD variables !
ALLOCATE(usstarL(imin:imax))
ALLOCATE(usstarR(imin:imax))
ALLOCATE(psstarL(imin:imax))
ALLOCATE(psstarR(imin:imax))

! NM !
ALLOCATE(epsL(-2:nx+3,-2:ny+3,-2:nz+3))
ALLOCATE(epsR(-2:nx+3,-2:ny+3,-2:nz+3))
ALLOCATE(csL(-2:nx+3,-2:ny+3,-2:nz+3))
ALLOCATE(csR(-2:nx+3,-2:ny+3,-2:nz+3))
ALLOCATE(fluxL(imin:imax,-2:nx+3,-2:ny+3,-2:nz+3))
ALLOCATE(fluxR(imin:imax,-2:nx+3,-2:ny+3,-2:nz+3))
ALLOCATE(uL(imin:imax,-2:nx+3,-2:ny+3,-2:nz+3))
ALLOCATE(uR(imin:imax,-2:nx+3,-2:ny+3,-2:nz+3))
ALLOCATE(primL(imin:imax,-2:nx+3,-2:ny+3,-2:nz+3))
ALLOCATE(primR(imin:imax,-2:nx+3,-2:ny+3,-2:nz+3))

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The LF flux, see for example. shu et al. for the original WENO paper !					    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LFNM(dir_in)
USE DEFINITION 
IMPLICIT NONE

! integer !
INTEGER, INTENT(IN) :: dir_in

! Integer !
INTEGER :: i, j, k, l

! Signal speeds !
REAL*8 :: cfsL, cfsR
REAL*8 :: sL, sR, splus

! Integer !
INTEGER :: ivn
INTEGER :: kx, ky, kz
INTEGER :: ibn, ibt1, ibt2

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
	ibn = ibx 
	ibt1 = iby
	ibt2 = ibz
	ivn = ivx
	kx = 0
ELSEIF(dir_in == y_dir) THEN
	ibn = iby
	ibt1 = ibz
	ibt2 = ibx
	ivn = ivy
	ky = 0
ELSEIF(dir_in == z_dir) THEN
	ibn = ibz
	ibt1 = ibx
	ibt2 = iby
	ivn = ivz
	kz = 0
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For NM !
!$OMP PARALLEL DO PRIVATE(cfsL, cfsR, sL, sR, splus) COLLAPSE(3) SCHEDULE(STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) PRIVATE(cfsL, cfsR, sL, sR, splus)
DO l = kz, nz 
	DO k = ky, ny
		DO j = kx, nx

			! Signal speed !
			cfsL = compute_signalspeed(csL(j,k,l), primL(ibn,j,k,l), primL(ibt1,j,k,l), primL(ibt2,j,k,l), primL(irho,j,k,l))
			cfsR = compute_signalspeed(csR(j,k,l), primR(ibn,j,k,l), primR(ibt1,j,k,l), primR(ibt2,j,k,l), primR(irho,j,k,l))
			sL = ABS(primL(ivn,j,k,l)) + cfsL
			sR = ABS(primR(ivn,j,k,l)) + cfsR
			splus = MAX(sL, sR)

			! fluxes !
			flux(imin:imax,j,k,l) = 0.5D0 * (fluxL (imin:imax,j,k,l) + fluxR (imin:imax,j,k,l) - splus * (uR(imin:imax,j,k,l) - uL (imin:imax,j,k,l)))

		END DO
	END DO
END DO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'lf = ', REAL(time_end - time_start) / rate
#endif

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The HLL flux, see for example, TORO's book on numerical hydrodynamics !					    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE HLLNM(dir_in)
USE DEFINITION 
IMPLICIT NONE

! integer !
INTEGER, INTENT(IN) :: dir_in

! Integer !
INTEGER :: i, j, k, l

! Signal speeds !
REAL*8 :: sL, sR
REAL*8 :: cfsL, cfsR
REAL*8 :: ubar, cbar

! Integer !
INTEGER :: ivn
INTEGER :: kx, ky, kz
INTEGER :: ibn, ibt1, ibt2

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
	ibn = ibx 
	ibt1 = iby
	ibt2 = ibz
	ivn = ivx
	kx = 0
ELSEIF(dir_in == y_dir) THEN
	ibn = iby
	ibt1 = ibz
	ibt2 = ibx
	ivn = ivy
	ky = 0
ELSEIF(dir_in == z_dir) THEN
	ibn = ibz
	ibt1 = ibx
	ibt2 = iby
	ivn = ivz
	kz = 0
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
! For NM !
!$OMP PARALLEL DO PRIVATE (sL, sR, cfsL, cfsR, ubar, cbar) COLLAPSE(3) SCHEDULE(STATIC) 
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) PRIVATE (sL, sR, cfsL, cfsR, ubar, cbar) 
DO l = kz, nz 
	DO k = ky, ny
		DO j = kx, nx

			! Signal speed !
			cfsL = compute_signalspeed(csL(j,k,l), primL(ibn,j,k,l), primL(ibt1,j,k,l), primL(ibt2,j,k,l), primL(irho,j,k,l))
			cfsR = compute_signalspeed(csR(j,k,l), primR(ibn,j,k,l), primR(ibt1,j,k,l), primR(ibt2,j,k,l), primR(irho,j,k,l))
			ubar = compute_roe(primL(ivn,j,k,l),primR(ivn,j,k,l),primL(irho,j,k,l),primR(irho,j,k,l))
			cbar = compute_roe(cfsL,cfsR,primL(irho,j,k,l),primR(irho,j,k,l))
			sL = min(primL(ivn,j,k,l) - cfsL, ubar - cbar)
			sR = max(primR(ivn,j,k,l) + cfsR, ubar + cbar)

			! Find the flux !
			IF(sL >= 0.0D0) THEN
				flux(imin:imax,j,k,l) = fluxL(imin:imax,j,k,l)
			ELSEIF(sL <= 0.0D0 .AND. sR >= 0.0D0) THEN
				DO i = imin, imax
					flux(i,j,k,l) = compute_fluxhll(fluxL(i,j,k,l),fluxR(i,j,k,l),uL(i,j,k,l),uR(i,j,k,l),sL,sR)
				END DO
			ELSEIF(sR <= 0.0D0) THEN
				flux(imin:imax,j,k,l) = fluxR(imin:imax,j,k,l)
			END IF

		END DO
	END DO
END DO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'hll = ', REAL(time_end - time_start) / rate
#endif

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The HLLC riemann solver, see Toro. 2009 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE HLLCNM(dir_in)
USE DEFINITION 
IMPLICIT NONE

! integer !
INTEGER, INTENT(IN) :: dir_in

! Integer !
INTEGER :: i, j, k, l

! Signal speeds !
REAL*8 :: sL, sR
REAL*8 :: pLs, pRs
REAL*8 :: pLt, pRt
REAL*8 :: b2L, b2R
REAL*8 :: vb, vbhll
REAL*8 :: cfsL, cfsR
REAL*8 :: ubar, cbar, sstar

! Integer !
INTEGER :: kx, ky, kz
INTEGER :: ivn, ivt1, ivt2
INTEGER :: ibn, ibt1, ibt2

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
	ibn = ibx 
	ibt1 = iby
	ibt2 = ibz
	ivn = ivx
	ivt1 = ivy
	ivt2 = ivz
	kx = 0
ELSEIF(dir_in == y_dir) THEN
	ibn = iby
	ibt1 = ibz
	ibt2 = ibx
	ivn = ivy
	ivt1 = ivz
	ivt2 = ivx
	ky = 0
ELSEIF(dir_in == z_dir) THEN
	ibn = ibz
	ibt1 = ibx
	ibt2 = iby
	ivn = ivz
	ivt1 = ivx
	ivt2 = ivy
	kz = 0
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For NM !
!$OMP PARALLEL DO PRIVATE (pLt, pRt, b2L, b2R, sL, sR, cfsL, cfsR, ubar, cbar, &
!$OMP sstar, pLs, pRs, ustarL, ustarR, u_hll, p_hll, vb, vbhll) COLLAPSE(3) SCHEDULE(STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) PRIVATE (pLt, pRt, b2L, b2R, &
!$ACC sL, sR, cfsL, cfsR, ubar, cbar, sstar, pLs, pRs, ustarL, ustarR, u_hll, p_hll, vb, vbhll)
DO l = kz, nz 
	DO k = ky, ny
		DO j = kx, nx

			! Signal speed !
			cfsL = compute_signalspeed(csL(j,k,l), primL(ibn,j,k,l), primL(ibt1,j,k,l), primL(ibt2,j,k,l), primL(irho,j,k,l))
			cfsR = compute_signalspeed(csR(j,k,l), primR(ibn,j,k,l), primR(ibt1,j,k,l), primR(ibt2,j,k,l), primR(irho,j,k,l))
			ubar = compute_roe(primL(ivn,j,k,l),primR(ivn,j,k,l),primL(irho,j,k,l),primR(irho,j,k,l))
			cbar = compute_roe(cfsL,cfsR,primL(irho,j,k,l),primR(irho,j,k,l))
			sL = min(primL(ivn,j,k,l) - cfsL, ubar - cbar)
			sR = max(primR(ivn,j,k,l) + cfsR, ubar + cbar)

			! Compute state accordingly !
			IF(sL >= 0.0d0) THEN
				flux(imin:imax,j,k,l) = fluxL(imin:imax,j,k,l)
			ELSEIF(sR <= 0.0d0) THEN
				flux(imin:imax,j,k,l) = fluxR(imin:imax,j,k,l)
			ELSE

				! bsquare !
				b2L = dot_product(primL(ibx:ibz,j,k,l), primL(ibx:ibz,j,k,l))
				b2R = dot_product(primR(ibx:ibz,j,k,l), primR(ibx:ibz,j,k,l))
				pLt = primL(itau,j,k,l) + 0.5D0*b2L
				pRt = primR(itau,j,k,l) + 0.5D0*b2R

				! Middle speed !
				sstar = compute_sstar_hllc(pLt,pRt,primL(irho,j,k,l),primR(irho,j,k,l),primL(ivn,j,k,l), &
																		primR(ivn,j,k,l),sL,sR,primL(ibn,j,k,l),primR(ibn,j,k,l))															

				! Get HLL state !
				u_hll(ibn) = primL(ibn,j,k,l) ! or primR, doesn't matter !
				u_hll(ibt1) = compute_hll(primL(ibt1,j,k,l),primR(ibt1,j,k,l),fluxL(ibt1,j,k,l),fluxR(ibt1,j,k,l),sL,sR)
				u_hll(ibt2) = compute_hll(primL(ibt2,j,k,l),primR(ibt2,j,k,l),fluxL(ibt2,j,k,l),fluxR(ibt2,j,k,l),sL,sR)
				u_hll(irho) = compute_hll(primL(irho,j,k,l),primR(irho,j,k,l),fluxL(irho,j,k,l),fluxR(irho,j,k,l),sL,sR)
				u_hll(ivt1) = compute_hll(uL(ivt1,j,k,l),uR(ivt1,j,k,l),fluxL(ivt1,j,k,l),fluxR(ivt1,j,k,l),sL,sR)
				u_hll(ivt2) = compute_hll(uL(ivt2,j,k,l),uR(ivt2,j,k,l),fluxL(ivt2,j,k,l),fluxR(ivt2,j,k,l),sL,sR)
				p_hll(ivn) = sstar
				p_hll(ivt1) = u_hll(ivt1)/u_hll(irho)
				p_hll(ivt2) = u_hll(ivt2)/u_hll(irho)

				! Star state !
				IF(sstar >= 0.0D0) THEN

					! initialize star state !
					ustarL(irho) = uL(irho,j,k,l)*(sL - primL(ivn,j,k,l))/(sL - sstar)
					ustarL(itau+1:ibx-1) = uL(itau+1:ibx-1,j,k,l)*(sL - primL(ivn,j,k,l))/(sL - sstar)

					! dot product !
					vb = dot_product(primL(ibx:ibz,j,k,l), primL(ivx:ivz,j,k,l))
					vbhll = dot_product(u_hll(ibx:ibz), p_hll(ivx:ivz))
		
					! Star pressure !
					pLs = primL(irho,j,k,l)*(sL - primL(ivn,j,k,l))*(sstar - primL(ivn,j,k,l)) + pLt
					
					! left star state, note we follow castro to treat rhooeps as a scalar !
					ustarL(ibn) = u_hll(ibn)
					ustarL(ibt1) = u_hll(ibt1)
					ustarL(ibt2) = u_hll(ibt2)
					ustarL(ivn) = ustarL(irho)*sstar
					ustarL(ivt1) = (uL(ivt1,j,k,l)*(sL - primL(ivn,j,k,l)) - (ustarL(ibn)*ustarL(ibt1) - primL(ibn,j,k,l)*primL(ibt1,j,k,l)))/(sL - sstar)
					ustarL(ivt2) = (uL(ivt2,j,k,l)*(sL - primL(ivn,j,k,l)) - (ustarL(ibn)*ustarL(ibt2) - primL(ibn,j,k,l)*primL(ibt2,j,k,l)))/(sL - sstar)
					ustarL(itau) = (uL(itau,j,k,l)*(sL - primL(ivn,j,k,l)) + (pLs*sstar - pLt*primL(ivn,j,k,l)) - ((vbhll*ustarL(ibn) - vb*primL(ibn,j,k,l))))/(sL - sstar)

					! Compute fluxes !
					flux(imin:imax,j,k,l) = fluxL(imin:imax,j,k,l) + sL*(ustarL(imin:imax) - uL(imin:imax,j,k,l))

				ELSE

					! initialize star state !
					ustarR(irho) = uR(irho,j,k,l)*(sR - primR(ivn,j,k,l))/(sR - sstar)
					ustarR(itau+1:ibx-1) = uR(itau+1:ibx-1,j,k,l)*(sR - primR(ivn,j,k,l))/(sR - sstar)

					! dot product !
					vb = dot_product(primR(ibx:ibz,j,k,l), primR(ivx:ivz,j,k,l))
					vbhll = dot_product(u_hll(ibx:ibz), p_hll(ivx:ivz))
		
					! Star pressure !
					pRs = primR(irho,j,k,l)*(sR - primR(ivn,j,k,l))*(sstar - primR(ivn,j,k,l)) + pRt
					
					! right star state, note we follow castro to treat rhooeps as a scalar !
					ustarR(ibn) = u_hll(ibn)
					ustarR(ibt1) = u_hll(ibt1)
					ustarR(ibt2) = u_hll(ibt2)
					ustarR(ivn) = ustarR(irho)*sstar
					ustarR(ivt1) = (uR(ivt1,j,k,l)*(sR - primR(ivn,j,k,l)) - (ustarR(ibn)*ustarR(ibt1) - primR(ibn,j,k,l)*primR(ibt1,j,k,l)))/(sR - sstar)
					ustarR(ivt2) = (uR(ivt2,j,k,l)*(sR - primR(ivn,j,k,l)) - (ustarR(ibn)*ustarR(ibt2) - primR(ibn,j,k,l)*primR(ibt2,j,k,l)))/(sR - sstar)
					ustarR(itau) = (uR(itau,j,k,l)*(sR - primR(ivn,j,k,l)) + (pRs*sstar - pRt*primR(ivn,j,k,l)) - ((vbhll*ustarR(ibn) - vb*primR(ibn,j,k,l))))/(sR - sstar)					

					! Compute fluxes !
					flux(imin:imax,j,k,l) = fluxR(imin:imax,j,k,l) + sR*(ustarR(imin:imax) - uR(imin:imax,j,k,l))

				END IF

			END IF

		END DO
	END DO
END DO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'hllc = ', REAL(time_end - time_start) / rate
#endif

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The HLLD riemann solver, 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE HLLDNM(dir_in)
USE DEFINITION 
IMPLICIT NONE

! integer !
INTEGER, INTENT(IN) :: dir_in

! Integer !
INTEGER :: i, j, k, l

! Signal speeds !
REAL*8 :: b2L, b2R
REAL*8 :: cfsl, cfsr
REAL*8 :: ubar, cbar
REAL*8 :: fac_1, fac_2
REAL*8 :: plt, prt, pls, prs
REAL*8 :: bdotu, bdotus, bdotuss
REAL*8 :: sl, sr, sm, sstarl, sstarr

! Integer !
INTEGER :: kx, ky, kz
INTEGER :: ivn, ivt1, ivt2
INTEGER :: ibn, ibt1, ibt2

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
	ibn = ibx 
	ibt1 = iby
	ibt2 = ibz
	ivn = ivx
	ivt1 = ivy
	ivt2 = ivz
	kx = 0
ELSEIF(dir_in == y_dir) THEN
	ibn = iby
	ibt1 = ibz
	ibt2 = ibx
	ivn = ivy
	ivt1 = ivz
	ivt2 = ivx
	ky = 0
ELSEIF(dir_in == z_dir) THEN
	ibn = ibz
	ibt1 = ibx
	ibt2 = iby
	ivn = ivz
	ivt1 = ivx
	ivt2 = ivy
	kz = 0
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For NM !
!$OMP PARALLEL DO PRIVATE (ubar, cbar, cfsl, cfsr, b2L, b2R, plt, prt, pls, prs, &
!$OMP sstarl, sstarr, bdotu, bdotus, bdotuss, sl, sr, sm, fac_1, fac_2, & 
!$OMP ustarL, ustarR, usstarL, usstarR, pstarL, pstarR, psstarL, psstarR) COLLAPSE(3) SCHEDULE(STATIC) 
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) PRIVATE (ubar, cbar, cfsl, cfsr, & 
!$ACC b2L, b2R, plt, prt, pls, prs, sstarl, sstarr, bdotu, bdotus, bdotuss, sl, sr, sm, fac_1, fac_2, & 
!$ACC ustarL, ustarR, usstarL, usstarR, pstarL, pstarR, psstarL, psstarR)
DO l = kz, nz 
	DO k = ky, ny
		DO j = kx, nx

			! Signal speed !
			cfsL = compute_signalspeed(csL(j,k,l), primL(ibn,j,k,l), primL(ibt1,j,k,l), primL(ibt2,j,k,l), primL(irho,j,k,l))
			cfsR = compute_signalspeed(csR(j,k,l), primR(ibn,j,k,l), primR(ibt1,j,k,l), primR(ibt2,j,k,l), primR(irho,j,k,l))
			ubar = compute_roe(primL(ivn,j,k,l),primR(ivn,j,k,l),primL(irho,j,k,l),primR(irho,j,k,l))
			cbar = compute_roe(cfsL,cfsR,primL(irho,j,k,l),primR(irho,j,k,l))
			sL = min(primL(ivn,j,k,l) - cfsL, ubar - cbar)
			sR = max(primR(ivn,j,k,l) + cfsR, ubar + cbar)

			! Compute state accordingly !
			IF(sL >= 0.0d0) THEN
				flux(imin:imax,j,k,l) = fluxL(imin:imax,j,k,l)
			ELSEIF(sR <= 0.0d0) THEN
				flux(imin:imax,j,k,l) = fluxR(imin:imax,j,k,l)
			ELSE

				! bsquare !
				b2L = dot_product(primL(ibx:ibz,j,k,l), primL(ibx:ibz,j,k,l))
				b2R = dot_product(primR(ibx:ibz,j,k,l), primR(ibx:ibz,j,k,l))
				pLt = primL(itau,j,k,l) + 0.5D0*b2L
				pRt = primR(itau,j,k,l) + 0.5D0*b2R

				! star-state wave speed !
				sM = compute_sstar_hllc(pLt,pRt,primL(irho,j,k,l),primR(irho,j,k,l),primL(ivn,j,k,l),&
																primR(ivn,j,k,l),sL,sR,primL(ibn,j,k,l),primR(ibn,j,k,l))

				! initialize star state !
				ustarL(irho:ibx-1) = uL(irho:ibx-1,j,k,l)*(sL - primL(ivn,j,k,l))/(sL - sM)
				ustarR(irho:ibx-1) = uR(irho:ibx-1,j,k,l)*(sR - primR(ivn,j,k,l))/(sR - sM)
				
				! Multi-state wave speed !
				sstarL = sM - ABS(primL(ibn,j,k,l))/DSQRT(ustarL(irho))
				sstarR = sM + ABS(primR(ibn,j,k,l))/DSQRT(ustarR(irho))

				! Condition to revert to HLLC, copied from pluto !
				IF (((sstarL - sL) <  1.0d-4*(sM - sL)) .OR. ((sstarR - sR) > - 1.0d-4*(sR - sM))) THEN

					! HLLC state !
					ustarL(ibn) = primL(ibn,j,k,l)
					ustarL(ibt1) = compute_hll(primL(ibt1,j,k,l),primR(ibt1,j,k,l),fluxL(ibt1,j,k,l),fluxR(ibt1,j,k,l),sL,sR)
					ustarL(ibt2) = compute_hll(primL(ibt2,j,k,l),primR(ibt2,j,k,l),fluxL(ibt2,j,k,l),fluxR(ibt2,j,k,l),sL,sR)
					ustarR(ibn) = primR(ibn,j,k,l)
					ustarR(ibt1) = ustarL(ibt1)
					ustarR(ibt2) = ustarL(ibt2)

					! Adjust speed !
					sstarL = sM 
					sstarR = SM

				ELSE

					! left states !
					fac_1 = (primL(irho,j,k,l)*(sL - primL(ivn,j,k,l))*(sL - sM) - primL(ibn,j,k,l)*primL(ibn,j,k,l))
					ustarL(ibn) = primL(ibn,j,k,l)
					ustarL(ibt1) = (primL(ibt1,j,k,l)*(primL(irho,j,k,l)*(sL - primL(ivn,j,k,l))*(sL - primL(ivn,j,k,l)) - primL(ibn,j,k,l)*primL(ibn,j,k,l)))/fac_1
					ustarL(ibt2) = (primL(ibt2,j,k,l)*(primL(irho,j,k,l)*(sL - primL(ivn,j,k,l))*(sL - primL(ivn,j,k,l)) - primL(ibn,j,k,l)*primL(ibn,j,k,l)))/fac_1

					! right states !
					fac_1 = (primR(irho,j,k,l)*(sR - primR(ivn,j,k,l))*(sR - sM) - primR(ibn,j,k,l)*primR(ibn,j,k,l))
					ustarR(ibn) = primR(ibn,j,k,l)
					ustarR(ibt1) = (primR(ibt1,j,k,l)*(primR(irho,j,k,l)*(sR - primR(ivn,j,k,l))*(sR - primR(ivn,j,k,l)) - primR(ibn,j,k,l)*primR(ibn,j,k,l)))/fac_1
					ustarR(ibt2) = (primR(ibt2,j,k,l)*(primR(irho,j,k,l)*(sR - primR(ivn,j,k,l))*(sR - primR(ivn,j,k,l)) - primR(ibn,j,k,l)*primR(ibn,j,k,l)))/fac_1

				END IF

				! Momentum !
				ustarL(ivn) = ustarL(irho)*sM
				ustarL(ivt1) = ustarL(ivt1) - (ustarL(ibn)*ustarL(ibt1) - primL(ibn,j,k,l)*primL(ibt1,j,k,l))/(sL - sM)
				ustarL(ivt2) = ustarL(ivt2) - (ustarL(ibn)*ustarL(ibt2) - primL(ibn,j,k,l)*primL(ibt2,j,k,l))/(sL - sM)

				! Momentum !
				ustarR(ivn) = ustarR(irho)*sM
				ustarR(ivt1) = ustarR(ivt1) - (ustarR(ibn)*ustarR(ibt1) - primR(ibn,j,k,l)*primR(ibt1,j,k,l))/(sR - sM)
				ustarR(ivt2) = ustarR(ivt2) - (ustarR(ibn)*ustarR(ibt2) - primR(ibn,j,k,l)*primR(ibt2,j,k,l))/(sR - sM)

				! Velocity !
				pstarL(ivn) = sM
				pstarL(ivt1) = ustarL(ivt1)/ustarL(irho)
				pstarL(ivt2) = ustarL(ivt2)/ustarL(irho)
				pstarR(ivn) = sM
				pstarR(ivt1) = ustarR(ivt1)/ustarR(irho)
				pstarR(ivt2) = ustarR(ivt2)/ustarR(irho)

				! Presure !
				pLs = primL(irho,j,k,l)*(sL - primL(ivn,j,k,l))*(sM - primL(ivn,j,k,l)) + pLt
				pRs = primR(irho,j,k,l)*(sR - primR(ivn,j,k,l))*(sM - primR(ivn,j,k,l)) + pRt
				
				! Internal energy !
				bdotu = dot_product(primL(ibx:ibz,j,k,l), primL(ivx:ivz,j,k,l))
				bdotus = dot_product(ustarL(ibx:ibz), pstarL(ivx:ivz))
				ustarL(itau) = ustarL(itau) + ((pLs*sM - pLt*primL(ivn,j,k,l)) - (bdotus - bdotu)*primL(ibn,j,k,l))/(sL - sM)

				! Right state !
				bdotu = dot_product(primR(ibx:ibz,j,k,l), primR(ivx:ivz,j,k,l))
				bdotus = dot_product(ustarR(ibx:ibz), pstarR(ivx:ivz))
				ustarR(itau) = ustarR(itau) + ((pRs*sM - pRt*primR(ivn,j,k,l)) - (bdotus - bdotu)*primR(ibn,j,k,l))/(sR - sM)

				! Choose state !
				IF(sstarL >= 0.0d0) THEN
					flux(imin:imax,j,k,l) = fluxL(imin:imax,j,k,l) + sL*(ustarL(imin:imax) - uL(imin:imax,j,k,l))
				ELSEIF(sstarR <= 0.0d0) THEN
					flux(imin:imax,j,k,l) = fluxR(imin:imax,j,k,l) + sR*(ustarR(imin:imax) - uR(imin:imax,j,k,l))
				ELSE

					! initialize star state !
					usstarL(irho) = ustarL(irho)
					usstarL(itau+1:ibx-1) = ustarL(itau+1:ibx-1)
					usstarR(irho) = ustarR(irho)
					usstarR(itau+1:ibx-1) = ustarR(itau+1:ibx-1)

					! factor !
					fac_1 = 1.0d0 / (DSQRT(ustarL(irho)) + DSQRT(ustarR(irho)))
					fac_2 = SIGN(1.0d0,primL(ibn,j,k,l))

					! left states !
					psstarL(ivn) = sM 
					psstarL(ivt1) = fac_1*(DSQRT(ustarL(irho))*pstarL(ivt1) + DSQRT(ustarR(irho))*pstarR(ivt1) + (ustarR(ibt1) - ustarL(ibt1))*fac_2)
					psstarL(ivt2) = fac_1*(DSQRT(ustarL(irho))*pstarL(ivt2) + DSQRT(ustarR(irho))*pstarR(ivt2) + (ustarR(ibt2) - ustarL(ibt2))*fac_2)
					usstarL(ibt1) = fac_1*(DSQRT(ustarL(irho))*ustarR(ibt1) + DSQRT(ustarR(irho))*ustarL(ibt1) + DSQRT(ustarL(irho)*ustarR(irho))*(pstarR(ivt1) - pstarL(ivt1))*fac_2)
					usstarL(ibt2) = fac_1*(DSQRT(ustarL(irho))*ustarR(ibt2) + DSQRT(ustarR(irho))*ustarL(ibt2) + DSQRT(ustarL(irho)*ustarR(irho))*(pstarR(ivt2) - pstarL(ivt2))*fac_2)

					! Momentum and internal energy !
					usstarL(ivn) = usstarL(irho)*psstarL(ivn)
					usstarL(ivt1) = usstarL(irho)*psstarL(ivt1)
					usstarL(ivt2) = usstarL(irho)*psstarL(ivt2)
					bdotus = dot_product(ustarL(ibx:ibz), pstarL(ivx:ivz))
					bdotuss = dot_product(usstarL(ibx:ibz), psstarL(ivx:ivz))
					usstarL(itau) = ustarL(itau) - DSQRT(ustarL(irho))*(bdotus - bdotuss)*fac_2
					
					! Right state
					psstarR(ivn) = sM
					psstarR(ivt1) = psstarL(ivt1)
					psstarR(ivt2) = psstarL(ivt2)
					usstarR(ibt1) = usstarL(ibt1)
					usstarR(ibt2) = usstarL(ibt2)

					! Momentum and internal energy !
					usstarR(ivn) = usstarR(irho)*psstarR(ivn)
					usstarR(ivt1) = usstarR(irho)*psstarR(ivt1)
					usstarR(ivt2) = usstarR(irho)*psstarR(ivt2)
					bdotus = dot_product(ustarR(ibx:ibz), pstarR(ivx:ivz))
					bdotuss = dot_product(usstarR(ibx:ibz), psstarR(ivx:ivz))
					usstarR(itau) = ustarR(itau) + DSQRT(ustarR(irho))*(bdotus - bdotuss)*fac_2

					! Choose state !
					IF(sM >= 0.0d0) THEN
						flux(imin:imax,j,k,l) = fluxL(imin:imax,j,k,l) + sL*(ustarL(imin:imax) - uL(imin:imax,j,k,l)) + sstarL*(usstarL(imin:imax) - ustarL(imin:imax))
					ELSE
						flux(imin:imax,j,k,l) = fluxR(imin:imax,j,k,l) + sR*(ustarR(imin:imax) - uR(imin:imax,j,k,l)) + sstarR*(usstarR(imin:imax) - ustarR(imin:imax))
					END IF

				END IF
			
			END IF

		END DO
	END DO
END DO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'hlld = ', REAL(time_end - time_start) / rate
#endif

END SUBROUTINE

END MODULE