!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This files contain all the riemann solvers available for !
! simulating hydrodynamics				   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE RIEMANN_MODULE
USE DEFINITION
IMPLICIT NONE

! Pressure !
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:,:) :: p1L, p1R

! Left and right hydro-states for DM/NM !
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:,:) :: eps1R, eps1L 
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:,:) :: eps2R, eps2L

! Speed of sound !
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:,:) :: cs1L, cs1R
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:,:) :: cs2L, cs2R

! Left and right fluxes, conserved quantity !
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: primL1, primR1
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: primL2, primR2
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: fluxL1, fluxR1, uL1, uR1
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: fluxL2, fluxR2, uL2, uR2

! the alpha in the LF flux !
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:) :: alpha1_x, alpha1_y, alpha1_z
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:) :: alpha2_x, alpha2_y, alpha2_z

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! function for computing the HLL flux
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	real(DP) function compute_fluxhll(yl,yr,xl,xr,sl,sr)
	implicit none
	real(DP) :: yl, yr, xl, xr, sl, sr
	compute_fluxhll = (sr*yl-sl*yr+sl*sr*(xr-xl))/(sr-sl)
	end function

	real(DP) function compute_roe(xl,xr,rhol,rhor)
	implicit none
	real(DP) :: xl,xr,rhol,rhor
	compute_roe = (sqrt(rhol)*xl + sqrt(rhor)*xr)/(sqrt(rhol) + sqrt(rhor))
	end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! function for computing the HLLC flux
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	real(DP) function compute_plr(pl,pr,rhol,rhor,vl,vr,sl,sr,star)
	implicit none
	real(DP) :: pl,pr,rhol,rhor,vl,vr,sl,sr,star
	compute_plr = 0.5D0*(pl+pr+rhol*(sl-vl)*(star-vl)+rhor*(sr-vr)*(star-vr))
	end function

	real(DP) function compute_sstar(pl,pr,rhol,rhor,vl,vr,sl,sr,vfl,vfr)
	implicit none
	real(DP) :: pl, pr, rhol, rhor, vl, vr, sl, sr, vfl, vfr
	compute_sstar = (pr - pl + rhol*vl*(sl-(vl-vfl)) - rhor*vr*(sr-(vr-vfr)))/&
			(rhol*(sl-(vl-vfl)) - rhor*(sr-(vr-vfr)))
	end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Assign left and right states and fluxes for riemann problem !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BUILDRIEMANN
USE DEFINITION 
IMPLICIT NONE

! Left and right fluxes, conserved quantity !
IF(DM_flag) THEN
	ALLOCATE(alpha1_x(imin1:imax1,-2:ny_1+3,-2:nz_1+3))
	ALLOCATE(alpha1_y(imin1:imax1,-2:ny_1+3,-2:nz_1+3))
	ALLOCATE(alpha1_z(imin1:imax1,-2:ny_1+3,-2:nz_1+3))
	ALLOCATE(p1L(1:n_dim,-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
	ALLOCATE(p1R(1:n_dim,-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
	ALLOCATE(eps1L(1:n_dim,-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
	ALLOCATE(eps1R(1:n_dim,-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
	ALLOCATE(cs1L(1:n_dim,-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
	ALLOCATE(cs1R(1:n_dim,-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
	ALLOCATE(fluxL1(imin1:imax1,1:n_dim,-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
	ALLOCATE(fluxR1(imin1:imax1,1:n_dim,-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
	ALLOCATE(uL1(imin1:imax1,1:n_dim,-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
	ALLOCATE(uR1(imin1:imax1,1:n_dim,-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
	ALLOCATE(primL1(imin1:imax1,1:n_dim,-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
	ALLOCATE(primR1(imin1:imax1,1:n_dim,-2:nx_1+3,-2:ny_1+3,-2:nz_1+3))
END IF

! NM !
ALLOCATE(alpha2_x(imin2:imax2,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(alpha2_y(imin2:imax2,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(alpha2_z(imin2:imax2,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(eps2L(1:n_dim,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(eps2R(1:n_dim,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(cs2L(1:n_dim,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(cs2R(1:n_dim,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(fluxL2(imin2:imax2,1:n_dim,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(fluxR2(imin2:imax2,1:n_dim,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(uL2(imin2:imax2,1:n_dim,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(uR2(imin2:imax2,1:n_dim,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(primL2(imin2:imax2,1:n_dim,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))
ALLOCATE(primR2(imin2:imax2,1:n_dim,-2:nx_2+3,-2:ny_2+3,-2:nz_2+3))

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The LF flux, see for example. shu et al. for the original WENO paper !					    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LFDM
USE DEFINITION 
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k, l, p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For DM !
DO CONCURRENT(j = nx_min_1 - 1:nx_part_1, k = ny_min_1 - 1:ny_part_1, l = nz_min_1 - 1:nz_part_1, i = imin1:imax1)
	flux_1(i,x_dir,j,k,l) = 0.5D0 * (fluxL1 (i,x_dir,j,k,l) + fluxR1 (i,x_dir,j,k,l) - alpha1_x(i,k,l) * (uR1(i,x_dir,j,k,l) - uL1 (i,x_dir,j,k,l)))
	IF(n_dim > 1) THEN
		flux_1(i,y_dir,j,k,l) = 0.5D0 * (fluxL1 (i,y_dir,j,k,l) + fluxR1 (i,y_dir,j,k,l) - alpha1_y(i,j,l) * (uR1 (i,y_dir,j,k,l) - uL1 (i,y_dir,j,k,l)))
	END IF
	IF(n_dim > 2) THEN
		flux_1(i,z_dir,j,k,l) = 0.5D0 * (fluxL1 (i,z_dir,j,k,l) + fluxR1 (i,z_dir,j,k,l) - alpha1_z(i,j,k) * (uR1 (i,z_dir,j,k,l) - uL1 (i,z_dir,j,k,l)))
	END IF
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The LF flux, see for example. shu et al. for the original WENO paper !					    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LFNM
USE DEFINITION 
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k, l, p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For DM !
DO CONCURRENT(j = nx_min_2 - 1:nx_part_2, k = ny_min_1 - 2:ny_part_2, l = nz_min_2 - 1:nz_part_2, i = imin2:imax2)
	flux_2(i,x_dir,j,k,l) = 0.5D0 * (fluxL2 (i,x_dir,j,k,l) + fluxR2 (i,x_dir,j,k,l) - alpha2_x(i,k,l) * (uR2(i,x_dir,j,k,l) - uL2 (i,x_dir,j,k,l)))
	IF(n_dim > 1) THEN
		flux_2(i,y_dir,j,k,l) = 0.5D0 * (fluxL2 (i,y_dir,j,k,l) + fluxR2 (i,y_dir,j,k,l) - alpha2_y(i,j,l) * (uR2 (i,y_dir,j,k,l) - uL2 (i,y_dir,j,k,l)))
	END IF
	IF(n_dim > 2) THEN
		flux_2(i,z_dir,j,k,l) = 0.5D0 * (fluxL2 (i,z_dir,j,k,l) + fluxR2 (i,z_dir,j,k,l) - alpha2_z(i,j,k) * (uR2 (i,z_dir,j,k,l) - uL2 (i,z_dir,j,k,l)))
	END IF
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

END MODULE
