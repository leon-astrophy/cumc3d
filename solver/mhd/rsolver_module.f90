!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This files contain all the riemann solvers available for !
! simulating hydrodynamics				   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE RIEMANN_MODULE
USE DEFINITION
IMPLICIT NONE

! HLLC solver variables  
REAL*8, ALLOCATABLE, DIMENSION(:) :: ustarL, ustarR
REAL*8, ALLOCATABLE, DIMENSION(:) :: fstarL, fstarR, dstar

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
	implicit none
	REAL*8 :: cs, bn, bt1, bt2, rho
	REAL*8 :: a2_mhd, bn_mhd, bt1_mhd, bt2_mhd, b2_mhd
	a2_mhd = cs**2
	bn_mhd = SQRT(bn**2/rho)
	bt1_mhd = SQRT(bt1**2/rho)
	bt2_mhd = SQRT(bt2**2/rho)
	b2_mhd = bn_mhd**2 + bt1_mhd**2 + bt2_mhd**2
	compute_signalspeed = SQRT(0.5D0*(a2_mhd + b2_mhd + SQRT((a2_mhd + b2_mhd)**2 - 4.0D0*a2_mhd*bn_mhd**2)))
	end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! function for computing the HLL flux
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	REAL*8 function compute_fluxhll(yl,yr,xl,xr,sl,sr)
	implicit none
	REAL*8 :: yl, yr, xl, xr, sl, sr
	compute_fluxhll = (sr*yl-sl*yr+sl*sr*(xr-xl))/(sr-sl)
	end function

	REAL*8 function compute_roe(xl,xr,rhol,rhor)
	implicit none
	REAL*8 :: xl,xr,rhol,rhor
	compute_roe = (sqrt(rhol)*xl + sqrt(rhor)*xr)/(sqrt(rhol) + sqrt(rhor))
	end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! function for computing the HLLC flux
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	REAL*8 function compute_plr(pl,pr,rhol,rhor,vl,vr,sl,sr,star)
	implicit none
	REAL*8 :: pl,pr,rhol,rhor,vl,vr,sl,sr,star
	compute_plr = 0.5D0*(pl+pr+rhol*(sl-vl)*(star-vl)+rhor*(sr-vr)*(star-vr))
	end function

	REAL*8 function compute_sstar(pl,pr,rhol,rhor,vl,vr,sl,sr)
	implicit none
	REAL*8 :: pl, pr, rhol, rhor, vl, vr, sl, sr
	compute_sstar = (pr - pl + rhol*vl*(sl-vl) - rhor*vr*(sr-vr))/&
			(rhol*(sl-vl) - rhor*(sr-vr))
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
ALLOCATE(fstarL(imin2:imax2))
ALLOCATE(fstarR(imin2:imax2))
ALLOCATE(dstar(imin2:imax2))

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
kx = 0
ky = 0
kz = 0

! Assign !
IF(dir_in == x_dir) THEN
	ibn = ibx 
	ibt1 = iby
	ibt2 = ibz
	ivn = ivel2_x
	kx = 1
ELSEIF(dir_in == y_dir) THEN
	ibn = iby
	ibt1 = ibz
	ibt2 = ibx
	ivn = ivel2_y
	ky = 1
ELSEIF(dir_in == z_dir) THEN
	ibn = ibz
	ibt1 = ibx
	ibt2 = iby
	ivn = ivel2_z
	kz = 1
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For NM !
!$OMP PARALLEL DO PRIVATE(cfsL, cfsR, sL, sR, splus) COLLAPSE(3) SCHEDULE(STATIC) 
DO l = nz_min_2 - kz, nz_part_2
	DO k = ny_min_2 - ky, ny_part_2
		DO j = nx_min_2 - kx, nx_part_2 

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign !
kx = 0
ky = 0
kz = 0

! Assign !
IF(dir_in == x_dir) THEN
	ibn = ibx 
	ibt1 = iby
	ibt2 = ibz
	ivn = ivel2_x
	kx = 1
ELSEIF(dir_in == y_dir) THEN
	ibn = iby
	ibt1 = ibz
	ibt2 = ibx
	ivn = ivel2_y
	ky = 1
ELSEIF(dir_in == z_dir) THEN
	ibn = ibz
	ibt1 = ibx
	ibt2 = iby
	ivn = ivel2_z
	kz = 1
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For NM !
!$OMP PARALLEL DO PRIVATE (sL, sR, cfsL, cfsR, ubar, cbar) COLLAPSE(3) SCHEDULE(STATIC) 
DO l = nz_min_2 - kz, nz_part_2
	DO k = ny_min_2 - ky, ny_part_2
		DO j = nx_min_2 - kx, nx_part_2

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
			ELSEIF(sL < 0.0D0 .AND. sR > 0.0D0) THEN
				DO i = imin2, imax2
					flux_2(i,j,k,l) = compute_fluxhll(fluxL2(i,j,k,l),fluxR2(i,j,k,l),uL2(i,j,k,l),uR2(i,j,k,l),sL,sR)
				END DO
			ELSEIF(sR <= 0.0D0) THEN
				flux_2(imin2:imax2,j,k,l) = fluxR2(imin2:imax2,j,k,l)
			END IF

		END DO
	END DO
END DO
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
REAL*8 :: cfsL, cfsR
REAL*8 :: ubar, cbar, sstar
REAL*8 :: pLs, pRs, pLR
REAL*8 :: pLt, pRt
REAL*8 :: b2L, b2R

! Integer !
INTEGER :: ivn
INTEGER :: kx, ky, kz
INTEGER :: ibn, ibt1, ibt2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign !
kx = 0
ky = 0
kz = 0

! Assign !
IF(dir_in == x_dir) THEN
	ibn = ibx 
	ibt1 = iby
	ibt2 = ibz
	ivn = ivel2_x
	kx = 1
ELSEIF(dir_in == y_dir) THEN
	ibn = iby
	ibt1 = ibz
	ibt2 = ibx
	ivn = ivel2_y
	ky = 1
ELSEIF(dir_in == z_dir) THEN
	ibn = ibz
	ibt1 = ibx
	ibt2 = iby
	ivn = ivel2_z
	kz = 1
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For NM !
!$OMP PARALLEL DO PRIVATE (pLt, pRt, b2L, b2R, sL, sR, cfsL, cfsR, ubar, cbar, sstar, pLs, pRs, pLR, ustarL, ustarR, fstarL, fstarR, dstar) COLLAPSE(3) SCHEDULE(STATIC) 
DO l = nz_min_2 - kz, nz_part_2
	DO k = ny_min_2 - ky, ny_part_2
		DO j = nx_min_2 - kx, nx_part_2

			! Signal speed !
			cfsL = compute_signalspeed(cs2L(j,k,l), primL2(ibn,j,k,l), primL2(ibt1,j,k,l), primL2(ibt2,j,k,l), primL2(irho2,j,k,l))
			cfsR = compute_signalspeed(cs2R(j,k,l), primR2(ibn,j,k,l), primR2(ibt1,j,k,l), primR2(ibt2,j,k,l), primR2(irho2,j,k,l))
			ubar = compute_roe(primL2(ivn,j,k,l),primR2(ivn,j,k,l),primL2(irho2,j,k,l),primR2(irho2,j,k,l))
			cbar = compute_roe(cfsL,cfsR,primL2(irho2,j,k,l),primR2(irho2,j,k,l))
			sL = min(primL2(ivn,j,k,l) - cfsL, ubar - cbar)
			sR = max(primR2(ivn,j,k,l) + cfsR, ubar + cbar)

			! bsquare !
			b2L = dot_product(primL2(ibx:ibz,j,k,l), primL2(ibx:ibz,j,k,l))
			b2R = dot_product(primR2(ibx:ibz,j,k,l), primR2(ibx:ibz,j,k,l))
			pLt = primL2(itau2,j,k,l) + 0.5D0*b2L
			pRt = primR2(itau2,j,k,l) + 0.5D0*b2R

			! Get wave speed in star region !
			sstar = compute_sstar(pLt, pRt, primL2(irho2,j,k,l), primR2(irho2,j,k,l), primL2(ivn,j,k,l), primR2(ivn,j,k,l), sL, sR)

			! star state !
			dstar(imin2:imax2) = 0.0D0
			dstar(ivn) = 1.0D0
			dstar(itau2) = sstar
			IF(dual_energy == 1) THEN
				dstar(ieps2) = sstar
			END IF

			! Get star pressure !
			pLs = pLt + primL2(irho2,j,k,l)*(sstar - primL2(ivn,j,k,l))*(sL - primL2(ivn,j,k,l))
			pRs = pRt + primR2(irho2,j,k,l)*(sstar - primR2(ivn,j,k,l))*(sR - primR2(ivn,j,k,l))

			! Assign star state of conservative variables !
			ustarL(imin2:imax2) = (sL*uL2(imin2:imax2,j,k,l) - fluxL2(imin2:imax2,j,k,l) + pLs*dstar(imin2:imax2))/(sL - sstar)
			ustarR(imin2:imax2) = (sR*uR2(imin2:imax2,j,k,l) - fluxR2(imin2:imax2,j,k,l) + pRs*dstar(imin2:imax2))/(sR - sstar)

			! Star state fluxes !
			fstarR(imin2:imax2) = fluxR2(imin2:imax2,j,k,l) + sR*(ustarR(imin2:imax2) - uR2(imin2:imax2,j,k,l))
			fstarL(imin2:imax2) = fluxL2(imin2:imax2,j,k,l) + sL*(ustarL(imin2:imax2) - uL2(imin2:imax2,j,k,l))
		
			! Get the output flux !
			IF(sL >= 0.0D0) THEN
				flux_2(imin2:imax2,j,k,l) = fluxL2(imin2:imax2,j,k,l)
			ELSEIF(sL <= 0.0D0 .AND. sstar >= 0.0D0) THEN
				flux_2(imin2:imax2,j,k,l) = fstarL(imin2:imax2)
			ELSEIF(sstar <= 0.0D0 .AND. sR >= 0.0D0) THEN
				flux_2(imin2:imax2,j,k,l) = fstarR(imin2:imax2)
			ELSEIF(sR <= 0.0D0) THEN
				flux_2(imin2:imax2,j,k,l) = fluxR2(imin2:imax2,j,k,l)
			END IF

		END DO
	END DO
END DO
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

END MODULE
