!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This files contain all the riemann solvers available for !
! simulating hydrodynamics				   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE RIEMANN_MODULE
USE DEFINITION
IMPLICIT NONE

! the alpha in the LF flux !
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:) :: alpha1_x, alpha1_y, alpha1_z
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:) :: alpha2_x, alpha2_y, alpha2_z

! DM moving grid !
REAL (DP), ALLOCATABLE, DIMENSION(:) :: vf1xR, vf1xL
REAL (DP), ALLOCATABLE, DIMENSION(:) :: vf1yR, vf1yL
REAL (DP), ALLOCATABLE, DIMENSION(:) :: vf1zR, vf1zL

! NM moving grid !
REAL (DP), ALLOCATABLE, DIMENSION(:) :: vf2xR, vf2xL
REAL (DP), ALLOCATABLE, DIMENSION(:) :: vf2yR, vf2yL
REAL (DP), ALLOCATABLE, DIMENSION(:) :: vf2zR, vf2zL

! Left and right hydro-states for DM/NM !
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:,:) :: eps1R, eps1L 
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:,:) :: eps2R, eps2L

! Pressure !
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:,:) :: p1L, p1R

! Speed of sound !
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:,:) :: cs1L, cs1R
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:,:) :: cs2L, cs2R

! Left and right fluxes, conserved quantity !
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: primL1, primR1
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: primL2, primR2
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: fluxL1, fluxR1, uL1, uR1
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: fluxL2, fluxR2, uL2, uR2

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
IF(RUNDM_flag) THEN
	ALLOCATE(alpha1_x(-2:ny_1+3,-2:nz_1+3,imin1:imax1))
	ALLOCATE(alpha1_y(-2:nx_1+3,-2:nz_1+3,imin1:imax1))
	ALLOCATE(alpha1_z(-2:nx_1+3,-2:ny_1+3,imin1:imax1))

	IF(movinggriddm_flag) THEN
		ALLOCATE (vf1xL(-2:nx_1+3))
		ALLOCATE (vf1xR(-2:nx_1+3))
		ALLOCATE (vf1yL(-2:ny_1+3))
		ALLOCATE (vf1yR(-2:ny_1+3))
		ALLOCATE (vf1zL(-2:nz_1+3))
		ALLOCATE (vf1zR(-2:nz_1+3))
	END IF

	ALLOCATE(p1L(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3,1:n_dim))
	ALLOCATE(p1R(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3,1:n_dim))
	ALLOCATE(eps1L(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3,1:n_dim))
	ALLOCATE(eps1R(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3,1:n_dim))
	ALLOCATE(cs1L(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3,1:n_dim))
	ALLOCATE(cs1R(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3,1:n_dim))
	ALLOCATE(fluxL1(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3,imin1:imax1,1:n_dim))
	ALLOCATE(fluxR1(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3,imin1:imax1,1:n_dim))
	ALLOCATE(uL1(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3,imin1:imax1,1:n_dim))
	ALLOCATE(uR1(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3,imin1:imax1,1:n_dim))
	ALLOCATE(primL1(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3,imin1:imax1,1:n_dim))
	ALLOCATE(primR1(-2:nx_1+3,-2:ny_1+3,-2:nz_1+3,imin1:imax1,1:n_dim))
END IF

! NM !
ALLOCATE(alpha2_x(-2:ny_2+3,-2:nz_2+3,imin2:imax2))
ALLOCATE(alpha2_y(-2:nx_2+3,-2:nz_2+3,imin2:imax2))
ALLOCATE(alpha2_z(-2:nx_2+3,-2:ny_2+3,imin2:imax2))

! Moving grid !
IF(movinggridnm_flag) THEN
	ALLOCATE (vf2xL(-2:nx_2+3))
	ALLOCATE (vf2xR(-2:nx_2+3))
	ALLOCATE (vf2yL(-2:ny_2+3))
	ALLOCATE (vf2yR(-2:ny_2+3))
	ALLOCATE (vf2zL(-2:nz_2+3))
	ALLOCATE (vf2zR(-2:nz_2+3))
END IF

ALLOCATE(eps2L(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3,1:n_dim))
ALLOCATE(eps2R(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3,1:n_dim))
ALLOCATE(cs2L(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3,1:n_dim))
ALLOCATE(cs2R(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3,1:n_dim))
ALLOCATE(fluxL2(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3,imin2:imax2,1:n_dim))
ALLOCATE(fluxR2(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3,imin2:imax2,1:n_dim))
ALLOCATE(uL2(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3,imin2:imax2,1:n_dim))
ALLOCATE(uR2(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3,imin2:imax2,1:n_dim))
ALLOCATE(primL2(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3,imin2:imax2,1:n_dim))
ALLOCATE(primR2(-2:nx_2+3,-2:ny_2+3,-2:nz_2+3,imin2:imax2,1:n_dim))

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
	flux_1(j,k,l,i,x_dir) = 0.5D0 * (fluxL1 (j,k,l,i,x_dir) + fluxR1 (j,k,l,i,x_dir) - alpha1_x(k,l,i) * (uR1 (j,k,l,i,x_dir) - uL1 (j,k,l,i,x_dir)))
	IF(n_dim > 1) THEN
		flux_1(j,k,l,i,y_dir) = 0.5D0 * (fluxL1 (j,k,l,i,y_dir) + fluxR1 (j,k,l,i,y_dir) - alpha1_y(j,l,i) * (uR1 (j,k,l,i,y_dir) - uL1 (j,k,l,i,y_dir)))
	END IF
	IF(n_dim > 2) THEN
		flux_1(j,k,l,i,z_dir) = 0.5D0 * (fluxL1 (j,k,l,i,z_dir) + fluxR1 (j,k,l,i,z_dir) - alpha1_z(j,k,i) * (uR1 (j,k,l,i,z_dir) - uL1 (j,k,l,i,z_dir)))
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
	flux_2(j,k,l,i,x_dir) = 0.5D0 * (fluxL2 (j,k,l,i,x_dir) + fluxR2 (j,k,l,i,x_dir) - alpha2_x(k,l,i) * (uR2 (j,k,l,i,x_dir) - uL2 (j,k,l,i,x_dir)))
	IF(n_dim > 1) THEN
		flux_2(j,k,l,i,y_dir) = 0.5D0 * (fluxL2 (j,k,l,i,y_dir) + fluxR2 (j,k,l,i,y_dir) - alpha2_y(j,l,i) * (uR2 (j,k,l,i,y_dir) - uL2 (j,k,l,i,y_dir)))
	END IF
	IF(n_dim > 2) THEN
		flux_2(j,k,l,i,z_dir) = 0.5D0 * (fluxL2 (j,k,l,i,z_dir) + fluxR2 (j,k,l,i,z_dir) - alpha2_z(j,k,i) * (uR2 (j,k,l,i,z_dir) - uL2 (j,k,l,i,z_dir)))
	END IF
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

END MODULE