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

! Pressure !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: p1L, p1R

! Left and right hydro-states for DM/NM !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: eps1R, eps1L 
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: eps2R, eps2L

! Speed of sound !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: cs1L, cs1R
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: cs2L, cs2R

! Left and right fluxes, conserved quantity !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: primL1, primR1
REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: primL2, primR2
REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: fluxL1, fluxR1, uL1, uR1
REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: fluxL2, fluxR2, uL2, uR2

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! function for computing the signal speed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	REAL*8 function compute_signalspeed(cs, bn, bt1, bt2, rho)
	!$ACC ROUTINE SEQ
	implicit none
	REAL*8 :: cs, bn, bt1, bt2, rho
	REAL*8 :: a2_mhd, bn_mhd, bt1_mhd, bt2_mhd, b2_mhd
	a2_mhd = cs**2
	bn_mhd = DSQRT(bn**2/rho)
	bt1_mhd = DSQRT(bt1**2/rho)
	bt2_mhd = DSQRT(bt2**2/rho)
	b2_mhd = bn_mhd**2 + bt1_mhd**2 + bt2_mhd**2
	compute_signalspeed = DSQRT(0.5D0*(a2_mhd + b2_mhd + DSQRT(MAX((a2_mhd + b2_mhd)**2 - 4.0D0*a2_mhd*bn_mhd**2, 0.0d0))))
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
	compute_sstar_hllc = (rhor*vr*(sr - vr) - rhol*vl*(sl-vl) + pl - pr + br**2 - bl**2) / &
									(rhor*(sr-vr) - rhol*(sl-vl))
	end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Assign left and right states and fluxes for riemann problem !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BUILDRIEMANN
USE DEFINITION 
IMPLICIT NONE

! HLLC variables !
ALLOCATE(ustarL(imin2:imax2))
ALLOCATE(ustarR(imin2:imax2))
ALLOCATE(pstarL(imin2:imax2))
ALLOCATE(pstarR(imin2:imax2))
ALLOCATE(u_hll(imin2:imax2))
ALLOCATE(p_hll(imin2:imax2))

! HLLD variables !
ALLOCATE(usstarL(imin2:imax2))
ALLOCATE(usstarR(imin2:imax2))
ALLOCATE(psstarL(imin2:imax2))
ALLOCATE(psstarR(imin2:imax2))

! NM !
ALLOCATE(eps2L(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(eps2R(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(cs2L(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(cs2R(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(fluxL2(imin2:imax2,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(fluxR2(imin2:imax2,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(uL2(imin2:imax2,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(uR2(imin2:imax2,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(primL2(imin2:imax2,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(primR2(imin2:imax2,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))

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
	ivn = ivel2_x
	kx = 0
ELSEIF(dir_in == y_dir) THEN
	ibn = iby
	ibt1 = ibz
	ibt2 = ibx
	ivn = ivel2_y
	ky = 0
ELSEIF(dir_in == z_dir) THEN
	ibn = ibz
	ibt1 = ibx
	ibt2 = iby
	ivn = ivel2_z
	kz = 0
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For NM !
!$OMP PARALLEL DO PRIVATE(cfsL, cfsR, sL, sR, splus) COLLAPSE(3) SCHEDULE(STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) PRIVATE(cfsL, cfsR, sL, sR, splus)
DO l = nz_min_2 - 1, nz_part_2 + kz
	DO k = ny_min_2 - 1, ny_part_2 + ky
		DO j = nx_min_2 - 1, nx_part_2 + kx

			! Signal speed !
			cfsL = compute_signalspeed(cs2L(j,k,l), primL2(ibn,j,k,l), primL2(ibt1,j,k,l), primL2(ibt2,j,k,l), primL2(irho2,j,k,l))
			cfsR = compute_signalspeed(cs2R(j,k,l), primR2(ibn,j,k,l), primR2(ibt1,j,k,l), primR2(ibt2,j,k,l), primR2(irho2,j,k,l))
			sL = ABS(primL2(ivn,j,k,l)) + cfsL
			sR = ABS(primR2(ivn,j,k,l)) + cfsR
			splus = MAX(sL, sR)

			! fluxes !
			flux_2(imin2:imax2,j,k,l) = 0.5D0 * (fluxL2 (imin2:imax2,j,k,l) + fluxR2 (imin2:imax2,j,k,l) - splus * (uR2(imin2:imax2,j,k,l) - uL2 (imin2:imax2,j,k,l)))

		END DO
	END DO
END DO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
	ibn = ibx 
	ibt1 = iby
	ibt2 = ibz
	ivn = ivel2_x
	kx = 0
ELSEIF(dir_in == y_dir) THEN
	ibn = iby
	ibt1 = ibz
	ibt2 = ibx
	ivn = ivel2_y
	ky = 0
ELSEIF(dir_in == z_dir) THEN
	ibn = ibz
	ibt1 = ibx
	ibt2 = iby
	ivn = ivel2_z
	kz = 0
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
! For NM !
!$OMP PARALLEL DO PRIVATE (sL, sR, cfsL, cfsR, ubar, cbar) COLLAPSE(3) SCHEDULE(STATIC) 
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) PRIVATE (sL, sR, cfsL, cfsR, ubar, cbar) 
DO l = nz_min_2 - 1, nz_part_2 + kz
	DO k = ny_min_2 - 1, ny_part_2 + ky
		DO j = nx_min_2 - 1, nx_part_2 + kx

			! Signal speed !
			cfsL = compute_signalspeed(cs2L(j,k,l), primL2(ibn,j,k,l), primL2(ibt1,j,k,l), primL2(ibt2,j,k,l), primL2(irho2,j,k,l))
			cfsR = compute_signalspeed(cs2R(j,k,l), primR2(ibn,j,k,l), primR2(ibt1,j,k,l), primR2(ibt2,j,k,l), primR2(irho2,j,k,l))
			ubar = compute_roe(primL2(ivn,j,k,l),primR2(ivn,j,k,l),primL2(irho2,j,k,l),primR2(irho2,j,k,l))
			cbar = compute_roe(cfsL,cfsR,primL2(irho2,j,k,l),primR2(irho2,j,k,l))
			sL = min(primL2(ivn,j,k,l) - cfsL, ubar - cbar)
			sR = max(primR2(ivn,j,k,l) + cfsR, ubar + cbar)

			! Find the flux !
			IF(sL >= 0.0D0) THEN
				flux_2(imin2:imax2,j,k,l) = fluxL2(imin2:imax2,j,k,l)
			ELSEIF(sL <= 0.0D0 .AND. sR >= 0.0D0) THEN
				DO i = imin2, imax2
					flux_2(i,j,k,l) = compute_fluxhll(fluxL2(i,j,k,l),fluxR2(i,j,k,l),uL2(i,j,k,l),uR2(i,j,k,l),sL,sR)
				END DO
			ELSEIF(sR <= 0.0D0) THEN
				flux_2(imin2:imax2,j,k,l) = fluxR2(imin2:imax2,j,k,l)
			END IF

		END DO
	END DO
END DO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL system_clock(time_end)
#ifdef DEBUG
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
	ivn = ivel2_x
	ivt1 = ivel2_y
	ivt2 = ivel2_z
	kx = 0
ELSEIF(dir_in == y_dir) THEN
	ibn = iby
	ibt1 = ibz
	ibt2 = ibx
	ivn = ivel2_y
	ivt1 = ivel2_z
	ivt2 = ivel2_x
	ky = 0
ELSEIF(dir_in == z_dir) THEN
	ibn = ibz
	ibt1 = ibx
	ibt2 = iby
	ivn = ivel2_z
	ivt1 = ivel2_x
	ivt2 = ivel2_y
	kz = 0
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For NM !
!$OMP PARALLEL DO PRIVATE (pLt, pRt, b2L, b2R, sL, sR, cfsL, cfsR, ubar, cbar, &
!$OMP sstar, pLs, pRs, ustarL, ustarR, u_hll, p_hll, vb, vbhll) COLLAPSE(3) SCHEDULE(STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) PRIVATE (pLt, pRt, b2L, b2R, &
!$ACC sL, sR, cfsL, cfsR, ubar, cbar, sstar, pLs, pRs, ustarL, ustarR, u_hll, p_hll, vb, vbhll)
DO l = nz_min_2 - 1, nz_part_2 + kz
	DO k = ny_min_2 - 1, ny_part_2 + ky
		DO j = nx_min_2 - 1, nx_part_2 + kx

			! Signal speed !
			cfsL = compute_signalspeed(cs2L(j,k,l), primL2(ibn,j,k,l), primL2(ibt1,j,k,l), primL2(ibt2,j,k,l), primL2(irho2,j,k,l))
			cfsR = compute_signalspeed(cs2R(j,k,l), primR2(ibn,j,k,l), primR2(ibt1,j,k,l), primR2(ibt2,j,k,l), primR2(irho2,j,k,l))
			ubar = compute_roe(primL2(ivn,j,k,l),primR2(ivn,j,k,l),primL2(irho2,j,k,l),primR2(irho2,j,k,l))
			cbar = compute_roe(cfsL,cfsR,primL2(irho2,j,k,l),primR2(irho2,j,k,l))
			sL = min(primL2(ivn,j,k,l) - cfsL, ubar - cbar)
			sR = max(primR2(ivn,j,k,l) + cfsR, ubar + cbar)

			! Compute state accordingly !
			IF(sL >= 0.0d0) THEN
				flux_2(imin2:imax2,j,k,l) = fluxL2(imin2:imax2,j,k,l)
			ELSEIF(sR <= 0.0d0) THEN
				flux_2(imin2:imax2,j,k,l) = fluxR2(imin2:imax2,j,k,l)
			ELSE

				! bsquare !
				b2L = dot_product(primL2(ibx:ibz,j,k,l), primL2(ibx:ibz,j,k,l))
				b2R = dot_product(primR2(ibx:ibz,j,k,l), primR2(ibx:ibz,j,k,l))
				pLt = primL2(itau2,j,k,l) + 0.5D0*b2L
				pRt = primR2(itau2,j,k,l) + 0.5D0*b2R

				! Middle speed !
				sstar = compute_sstar_hllc(pLt,pRt,primL2(irho2,j,k,l),primR2(irho2,j,k,l),primL2(ivn,j,k,l), &
																		primR2(ivn,j,k,l),sL,sR,primL2(ibn,j,k,l),primR2(ibn,j,k,l))															

				! Get HLL state !
				u_hll(ibn) = primL2(ibn,j,k,l) ! or primR2, doesn't matter !
				u_hll(ibt1) = compute_hll(primL2(ibt1,j,k,l),primR2(ibt1,j,k,l),fluxL2(ibt1,j,k,l),fluxR2(ibt1,j,k,l),sL,sR)
				u_hll(ibt2) = compute_hll(primL2(ibt2,j,k,l),primR2(ibt2,j,k,l),fluxL2(ibt2,j,k,l),fluxR2(ibt2,j,k,l),sL,sR)
				u_hll(irho2) = compute_hll(primL2(irho2,j,k,l),primR2(irho2,j,k,l),fluxL2(irho2,j,k,l),fluxR2(irho2,j,k,l),sL,sR)
				u_hll(ivt1) = compute_hll(uL2(ivt1,j,k,l),uR2(ivt1,j,k,l),fluxL2(ivt1,j,k,l),fluxR2(ivt1,j,k,l),sL,sR)
				u_hll(ivt2) = compute_hll(uL2(ivt2,j,k,l),uR2(ivt2,j,k,l),fluxL2(ivt2,j,k,l),fluxR2(ivt2,j,k,l),sL,sR)
				p_hll(ivn) = sstar
				p_hll(ivt1) = u_hll(ivt1)/u_hll(irho2)
				p_hll(ivt2) = u_hll(ivt2)/u_hll(irho2)

				! Star state !
				IF(sstar >= 0.0D0) THEN

					! initialize star state !
					ustarL(imin2:imax2) = uL2(imin2:imax2,j,k,l)*(sL - primL2(ivn,j,k,l))/(sL - sstar)

					! dot product !
					vb = dot_product(primL2(ibx:ibz,j,k,l), primL2(ivel2_x:ivel2_z,j,k,l))
					vbhll = dot_product(u_hll(ibx:ibz), p_hll(ivel2_x:ivel2_z))
		
					! Star pressure !
					pLs = primL2(irho2,j,k,l)*(sL - primL2(ivn,j,k,l))*(sstar - primL2(ivn,j,k,l)) + pLt
					
					! left star state, note we follow castro to treat rhooeps as a scalar !
					ustarL(ibn) = u_hll(ibn)
					ustarL(ibt1) = u_hll(ibt1)
					ustarL(ibt2) = u_hll(ibt2)
					ustarL(ivn) = ustarL(irho2)*sstar
					ustarL(ivt1) = (uL2(ivt1,j,k,l)*(sL - primL2(ivn,j,k,l)) - (ustarL(ibn)*ustarL(ibt1) - primL2(ibn,j,k,l)*primL2(ibt1,j,k,l)))/(sL - sstar)
					ustarL(ivt2) = (uL2(ivt2,j,k,l)*(sL - primL2(ivn,j,k,l)) - (ustarL(ibn)*ustarL(ibt2) - primL2(ibn,j,k,l)*primL2(ibt2,j,k,l)))/(sL - sstar)
					ustarL(itau2) = (uL2(itau2,j,k,l)*(sL - primL2(ivn,j,k,l)) + (pLs*sstar - pLt*primL2(ivn,j,k,l)) - ((vbhll*ustarL(ibn) - vb*primL2(ibn,j,k,l))))/(sL - sstar)

					! Compute fluxes !
					flux_2(imin2:imax2,j,k,l) = fluxL2(imin2:imax2,j,k,l) + sL*(ustarL(imin2:imax2) - uL2(imin2:imax2,j,k,l))

				ELSE

					! initialize star state !
					ustarR(imin2:imax2) = uR2(imin2:imax2,j,k,l)*(sR - primR2(ivn,j,k,l))/(sR - sstar)

					! dot product !
					vb = dot_product(primR2(ibx:ibz,j,k,l), primR2(ivel2_x:ivel2_z,j,k,l))
					vbhll = dot_product(u_hll(ibx:ibz), p_hll(ivel2_x:ivel2_z))
		
					! Star pressure !
					pRs = primR2(irho2,j,k,l)*(sR - primR2(ivn,j,k,l))*(sstar - primR2(ivn,j,k,l)) + pRt
					
					! right star state, note we follow castro to treat rhooeps as a scalar !
					ustarR(ibn) = u_hll(ibn)
					ustarR(ibt1) = u_hll(ibt1)
					ustarR(ibt2) = u_hll(ibt2)
					ustarR(ivn) = ustarR(irho2)*sstar
					ustarR(ivt1) = (uR2(ivt1,j,k,l)*(sR - primR2(ivn,j,k,l)) - (ustarR(ibn)*ustarR(ibt1) - primR2(ibn,j,k,l)*primR2(ibt1,j,k,l)))/(sR - sstar)
					ustarR(ivt2) = (uR2(ivt2,j,k,l)*(sR - primR2(ivn,j,k,l)) - (ustarR(ibn)*ustarR(ibt2) - primR2(ibn,j,k,l)*primR2(ibt2,j,k,l)))/(sR - sstar)
					ustarR(itau2) = (uR2(itau2,j,k,l)*(sR - primR2(ivn,j,k,l)) + (pRs*sstar - pRt*primR2(ivn,j,k,l)) - ((vbhll*ustarR(ibn) - vb*primR2(ibn,j,k,l))))/(sR - sstar)					

					! Compute fluxes !
					flux_2(imin2:imax2,j,k,l) = fluxR2(imin2:imax2,j,k,l) + sR*(ustarR(imin2:imax2) - uR2(imin2:imax2,j,k,l))

				END IF

			END IF

		END DO
	END DO
END DO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
REAL*8 :: b2l, b2r
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
	ivn = ivel2_x
	ivt1 = ivel2_y
	ivt2 = ivel2_z
	kx = 0
ELSEIF(dir_in == y_dir) THEN
	ibn = iby
	ibt1 = ibz
	ibt2 = ibx
	ivn = ivel2_y
	ivt1 = ivel2_z
	ivt2 = ivel2_x
	ky = 0
ELSEIF(dir_in == z_dir) THEN
	ibn = ibz
	ibt1 = ibx
	ibt2 = iby
	ivn = ivel2_z
	ivt1 = ivel2_x
	ivt2 = ivel2_y
	kz = 0
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For NM !
!$OMP PARALLEL DO PRIVATE (ubar, cbar, cfsl, cfsr, b2l, b2r, plt, prt, pls, prs, &
!$OMP sstarl, sstarr, bdotu, bdotus, bdotuss, sl, sr, sm, fac_1, fac_2, & 
!$OMP ustarL, ustarR, usstarL, usstarR, pstarL, pstarR, psstarL, psstarR) COLLAPSE(3) SCHEDULE(STATIC) 
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) PRIVATE (ubar, cbar, cfsl, cfsr, & 
!$ACC b2l, b2r, plt, prt, pls, prs, sstarl, sstarr, bdotu, bdotus, bdotuss, sl, sr, sm, fac_1, fac_2, & 
!$ACC ustarL, ustarR, usstarL, usstarR, pstarL, pstarR, psstarL, psstarR)
DO l = nz_min_2 - 1, nz_part_2 + kz
	DO k = ny_min_2 - 1, ny_part_2 + ky
		DO j = nx_min_2 - 1, nx_part_2 + kx

			! Signal speed !
			cfsL = compute_signalspeed(cs2L(j,k,l), primL2(ibn,j,k,l), primL2(ibt1,j,k,l), primL2(ibt2,j,k,l), primL2(irho2,j,k,l))
			cfsR = compute_signalspeed(cs2R(j,k,l), primR2(ibn,j,k,l), primR2(ibt1,j,k,l), primR2(ibt2,j,k,l), primR2(irho2,j,k,l))
			ubar = compute_roe(primL2(ivn,j,k,l),primR2(ivn,j,k,l),primL2(irho2,j,k,l),primR2(irho2,j,k,l))
			cbar = compute_roe(cfsL,cfsR,primL2(irho2,j,k,l),primR2(irho2,j,k,l))
			sL = min(primL2(ivn,j,k,l) - cfsL, ubar - cbar)
			sR = max(primR2(ivn,j,k,l) + cfsR, ubar + cbar)

			! Compute state accordingly !
			IF(sL > 0.0d0) THEN
				flux_2(imin2:imax2,j,k,l) = fluxL2(imin2:imax2,j,k,l)
			ELSEIF(sR <= 0.0d0) THEN
				flux_2(imin2:imax2,j,k,l) = fluxR2(imin2:imax2,j,k,l)
			ELSE

				! bsquare !
				b2L = dot_product(primL2(ibx:ibz,j,k,l), primL2(ibx:ibz,j,k,l))
				b2R = dot_product(primR2(ibx:ibz,j,k,l), primR2(ibx:ibz,j,k,l))
				pLt = primL2(itau2,j,k,l) + 0.5D0*b2L
				pRt = primR2(itau2,j,k,l) + 0.5D0*b2R

				! star-state wave speed !
				sM = compute_sstar_hllc(pLt,pRt,primL2(irho2,j,k,l),primR2(irho2,j,k,l),primL2(ivn,j,k,l),&
																primR2(ivn,j,k,l),sL,sR,primL2(ibn,j,k,l),primR2(ibn,j,k,l))

				! initialize star state !
				ustarL(imin2:imax2) = uL2(imin2:imax2,j,k,l)*(sL - primL2(ivn,j,k,l))/(sL - sM)
				ustarR(imin2:imax2) = uR2(imin2:imax2,j,k,l)*(sR - primR2(ivn,j,k,l))/(sR - sM)
				
				! Multi-state wave speed !
				sstarL = sM - ABS(primL2(ibn,j,k,l))/DSQRT(ustarL(irho2))
				sstarR = sM + ABS(primR2(ibn,j,k,l))/DSQRT(ustarR(irho2))

				! Condition to revert to HLLC, copied from pluto !
				IF (((sstarL - sL) <  1.e-4*(sM - sL)) .OR. ((sstarR - sR) > -1.e-4*(sR - sM))) THEN

					! HLLC state !
					ustarL(ibn) = primL2(ibn,j,k,l)
					ustarL(ibt1) = compute_hll(primL2(ibt1,j,k,l),primR2(ibt1,j,k,l),fluxL2(ibt1,j,k,l),fluxR2(ibt1,j,k,l),sL,sR)
					ustarL(ibt2) = compute_hll(primL2(ibt2,j,k,l),primR2(ibt2,j,k,l),fluxL2(ibt2,j,k,l),fluxR2(ibt2,j,k,l),sL,sR)
					ustarR(ibn) = primR2(ibn,j,k,l)
					ustarR(ibt1) = ustarL(ibt1)
					ustarR(ibt2) = ustarL(ibt2)

					! Adjust speed !
					sstarL = sM 
					sstarR = SM

				ELSE

					! left states !
					fac_1 = (primL2(irho2,j,k,l)*(sL - primL2(ivn,j,k,l))*(sL - sM) - primL2(ibn,j,k,l)**2)
					ustarL(ibn) = primL2(ibn,j,k,l)
					ustarL(ibt1) = (primL2(ibt1,j,k,l)*(primL2(irho2,j,k,l)*(sL - primL2(ivn,j,k,l))**2 - primL2(ibn,j,k,l)**2))/fac_1
					ustarL(ibt2) = (primL2(ibt2,j,k,l)*(primL2(irho2,j,k,l)*(sL - primL2(ivn,j,k,l))**2 - primL2(ibn,j,k,l)**2))/fac_1

					! right states !
					fac_1 = (primR2(irho2,j,k,l)*(sR - primR2(ivn,j,k,l))*(sR - sM) - primR2(ibn,j,k,l)**2)
					ustarR(ibn) = primR2(ibn,j,k,l)
					ustarR(ibt1) = (primR2(ibt1,j,k,l)*(primR2(irho2,j,k,l)*(sR - primR2(ivn,j,k,l))**2 - primR2(ibn,j,k,l)**2))/fac_1
					ustarR(ibt2) = (primR2(ibt2,j,k,l)*(primR2(irho2,j,k,l)*(sR - primR2(ivn,j,k,l))**2 - primR2(ibn,j,k,l)**2))/fac_1

				END IF

				! Momentum !
				ustarL(ivn) = ustarL(irho2)*sM
				ustarL(ivt1) = ustarL(ivt1) - (ustarL(ibn)*ustarL(ibt1) - primL2(ibn,j,k,l)*primL2(ibt1,j,k,l))/(sL - sM)
				ustarL(ivt2) = ustarL(ivt2) - (ustarL(ibn)*ustarL(ibt2) - primL2(ibn,j,k,l)*primL2(ibt2,j,k,l))/(sL - sM)

				! Momentum !
				ustarR(ivn) = ustarR(irho2)*sM
				ustarR(ivt1) = ustarR(ivt1) - (ustarR(ibn)*ustarR(ibt1) - primR2(ibn,j,k,l)*primR2(ibt1,j,k,l))/(sR - sM)
				ustarR(ivt2) = ustarR(ivt2) - (ustarR(ibn)*ustarR(ibt2) - primR2(ibn,j,k,l)*primR2(ibt2,j,k,l))/(sR - sM)

				! Velocity !
				pstarL(ivn) = sM
				pstarL(ivt1) = ustarL(ivt1)/ustarL(irho2)
				pstarL(ivt2) = ustarL(ivt2)/ustarL(irho2)
				pstarR(ivn) = sM
				pstarR(ivt1) = ustarR(ivt1)/ustarR(irho2)
				pstarR(ivt2) = ustarR(ivt2)/ustarR(irho2)

				! Presure !
				pLs = primL2(irho2,j,k,l)*(sL - primL2(ivn,j,k,l))*(sM - primL2(ivn,j,k,l)) + pLt
				pRs = primR2(irho2,j,k,l)*(sR - primR2(ivn,j,k,l))*(sM - primR2(ivn,j,k,l)) + pRt
				
				! Internal energy !
				bdotu = dot_product(primL2(ibx:ibz,j,k,l), primL2(ivel2_x:ivel2_z,j,k,l))
				bdotus = dot_product(ustarL(ibx:ibz), pstarL(ivel2_x:ivel2_z))
				ustarL(itau2) = ustarL(itau2) + ((pLs*sM - pLt*primL2(ivn,j,k,l)) - (bdotus - bdotu)*primL2(ibn,j,k,l))/(sL - sM)

				! Right state !
				bdotu = dot_product(primR2(ibx:ibz,j,k,l), primR2(ivel2_x:ivel2_z,j,k,l))
				bdotus = dot_product(ustarR(ibx:ibz), pstarR(ivel2_x:ivel2_z))
				ustarR(itau2) = ustarR(itau2) + ((pRs*sM - pRt*primR2(ivn,j,k,l)) - (bdotus - bdotu)*primR2(ibn,j,k,l))/(sR - sM)

				! Choose state !
				IF(sstarL > 0.0d0) THEN
					flux_2(imin2:imax2,j,k,l) = fluxL2(imin2:imax2,j,k,l) + sL*(ustarL(imin2:imax2) - uL2(imin2:imax2,j,k,l))
				ELSEIF(sstarR <= 0.0d0) THEN
					flux_2(imin2:imax2,j,k,l) = fluxR2(imin2:imax2,j,k,l) + sR*(ustarR(imin2:imax2) - uR2(imin2:imax2,j,k,l))
				ELSE

					! initialize star state !
					usstarL(imin2:imax2) = ustarL(imin2:imax2)
					usstarR(imin2:imax2) = ustarR(imin2:imax2)

					! factor !
					fac_1 = 1.0d0 / (DSQRT(ustarL(irho2)) + DSQRT(ustarR(irho2)))
					fac_2 = SIGN(1.0d0,primL2(ibn,j,k,l))

					! left states !
					psstarL(ivn) = sM 
					psstarL(ivt1) = fac_1*(DSQRT(ustarL(irho2))*pstarL(ivt1) + DSQRT(ustarR(irho2))*pstarR(ivt1) + (ustarR(ibt1) - ustarL(ibt1))*fac_2)
					psstarL(ivt2) = fac_1*(DSQRT(ustarL(irho2))*pstarL(ivt2) + DSQRT(ustarR(irho2))*pstarR(ivt2) + (ustarR(ibt2) - ustarL(ibt2))*fac_2)
					usstarL(ibt1) = fac_1*(DSQRT(ustarL(irho2))*ustarR(ibt1) + DSQRT(ustarR(irho2))*ustarL(ibt1) + DSQRT(ustarL(irho2)*ustarR(irho2))*(pstarR(ivt1) - pstarL(ivt1))*fac_2)
					usstarL(ibt2) = fac_1*(DSQRT(ustarL(irho2))*ustarR(ibt2) + DSQRT(ustarR(irho2))*ustarL(ibt2) + DSQRT(ustarL(irho2)*ustarR(irho2))*(pstarR(ivt2) - pstarL(ivt2))*fac_2)

					! Momentum and internal energy !
					usstarL(ivn) = usstarL(irho2)*psstarL(ivn)
					usstarL(ivt1) = usstarL(irho2)*psstarL(ivt1)
					usstarL(ivt2) = usstarL(irho2)*psstarL(ivt2)
					bdotus = dot_product(ustarL(ibx:ibz), pstarL(ivel2_x:ivel2_z))
					bdotuss = dot_product(usstarL(ibx:ibz), psstarL(ivel2_x:ivel2_z))
					usstarL(itau2) = ustarL(itau2) - DSQRT(ustarL(irho2))*(bdotus - bdotuss)*fac_2
					
					! Right state
					psstarR(ivn) = sM
					psstarR(ivt1) = psstarL(ivt1)
					psstarR(ivt2) = psstarL(ivt2)
					usstarR(ibt1) = usstarL(ibt1)
					usstarR(ibt2) = usstarL(ibt2)

					! Momentum and internal energy !
					usstarR(ivn) = usstarR(irho2)*psstarR(ivn)
					usstarR(ivt1) = usstarR(irho2)*psstarR(ivt1)
					usstarR(ivt2) = usstarR(irho2)*psstarR(ivt2)
					bdotus = dot_product(ustarR(ibx:ibz), pstarR(ivel2_x:ivel2_z))
					bdotuss = dot_product(usstarR(ibx:ibz), psstarR(ivel2_x:ivel2_z))
					usstarR(itau2) = ustarR(itau2) + DSQRT(ustarR(irho2))*(bdotus - bdotuss)*fac_2

					! Choose state !
					IF(sM > 0.0d0) THEN
						flux_2(imin2:imax2,j,k,l) = fluxL2(imin2:imax2,j,k,l) + sL*(ustarL(imin2:imax2) - uL2(imin2:imax2,j,k,l)) + sstarL*(usstarL(imin2:imax2) - ustarL(imin2:imax2))
					ELSE
						flux_2(imin2:imax2,j,k,l) = fluxR2(imin2:imax2,j,k,l) + sR*(ustarR(imin2:imax2) - uR2(imin2:imax2,j,k,l)) + sstarR*(usstarR(imin2:imax2) - ustarR(imin2:imax2))
					END IF

				END IF
			
			END IF

		END DO
	END DO
END DO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

END MODULE