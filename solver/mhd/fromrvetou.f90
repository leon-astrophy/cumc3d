!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine converts the primitive variables to
! conservative variables (or vice versa)
! Written by Leung Shing Chi in 2016
! If you add your own variables in the WENO scheme,
! add your conversion step here.
! This subroutine takes in the U array and conversion mode
! Mode 0: From primitive to conservative
! Mode 1: From conservative to primitive
!
! Here is a reminder in how to add new physics:
! 1. Add your own module that contains the physicsb
! 2. Remind BuildWENO to include your quantity
! 3. Add the conversion here
! 4. Write a section in how to calculate the flux in Spatial
! 5. Add a flag to give signal to the program whenever you use the code
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE FROMRVETOU
USE DEFINITION
USE MHD_MODULE
IMPLICIT NONE

! Dummy variables
INTEGER :: i, j, k, l

! Squared velocity !
REAL*8 :: vsquare, bsquare
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NM Sector !

! Convert NM hydro
!$OMP PARALLEL PRIVATE(vsquare, bsquare)
!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
DO l = nz_min_2, nz_part_2
  DO k = ny_min_2, ny_part_2
    DO j = nx_min_2, nx_part_2
      vsquare = dot_product(prim2(ivel2_x:ivel2_z,j,k,l), prim2(ivel2_x:ivel2_z,j,k,l))
      bsquare = dot_product(prim2(ibx:ibz,j,k,l), prim2(ibx:ibz,j,k,l))
		  cons2(imin2:imax2,j,k,l) = prim2(imin2:imax2,j,k,l)*prim2(irho2,j,k,l)
	    cons2(irho2,j,k,l) = prim2(irho2,j,k,l)
	    cons2(itau2,j,k,l) = prim2(irho2,j,k,l)*(epsilon2(j,k,l) + 0.5D0*vsquare) + 0.5D0*bsquare
		  cons2(ibx:ibz,j,k,l) = prim2(ibx:ibz,j,k,l)
	  END DO
  END DO
END DO
!$OMP END DO

! Dual energy !
IF(dual_energy) THEN
  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  DO l = nz_min_2, nz_part_2
    DO k = ny_min_2, ny_part_2
      DO j = nx_min_2, nx_part_2
		    cons2(ieps2,j,k,l) = prim2(irho2,j,k,l) * epsilon2(j,k,l)
	    END DO
	  END DO
  END DO
  !$OMP END DO
END IF
!$OMP END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Convert back to primitive variables !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FROMUTORVE
USE DEFINITION
USE MHD_MODULE
IMPLICIT NONE

! Dummy variables
INTEGER :: i, j, k, l

! Dummy variables !
REAL*8 :: et_2, bige_2

! Squared velocity !
REAL*8 :: vsquare, bsquare

! For epsilon !
REAL*8 :: diff, factor, dummy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For NM sectors !

! Convert the NM hydro
!$OMP PARALLEL PRIVATE(vsquare, bsquare, diff, factor, et_2, bige_2, dummy)
!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
DO l = nz_min_2, nz_part_2
  DO k = ny_min_2, ny_part_2
    DO j = nx_min_2, nx_part_2
		  prim2(imin2:imax2,j,k,l) = cons2(imin2:imax2,j,k,l)/cons2(irho2,j,k,l)
	    prim2(irho2,j,k,l) = cons2(irho2,j,k,l)
      prim2(ibx:ibz,j,k,l) = cons2(ibx:ibz,j,k,l)
      bsquare = dot_product(prim2(ibx:ibz,j,k,l), prim2(ibx:ibz,j,k,l))
      cons2(itau2,j,k,l) = cons2(itau2,j,k,l) - 0.5D0*bsquare
	  END DO
  END DO
END DO
!$OMP END DO

! Determine the epsilon for epsilon equation !
IF (dual_energy) THEN
  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  DO l = nz_min_2, nz_part_2
    DO k = ny_min_2, ny_part_2
      DO j = nx_min_2, nx_part_2
		    prim2(ieps2,j,k,l) = cons2(ieps2,j,k,l)
		    bige_2 = cons2(itau2,j,k,l)/cons2(irho2,j,k,l)
        vsquare = dot_product(prim2(ivel2_x:ivel2_z,j,k,l), prim2(ivel2_x:ivel2_z,j,k,l))
		    et_2 = bige_2 - 0.5D0 * vsquare
        diff = et_2 - 1.0D-4*bige_2
        factor = MAX(SIGN(1.0D0, diff), 0.0D0)
			  epsilon2(j,k,l) = factor*et_2 + (1.0d0 - factor)*cons2 (ieps2,j,k,l)/cons2(irho2,j,k,l)
	    END DO
    END DO
  END DO
  !$OMP END DO
ELSE
  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  DO l = nz_min_2, nz_part_2
    DO k = ny_min_2, ny_part_2
      DO j = nx_min_2, nx_part_2
        bige_2 = cons2(itau2,j,k,l)/cons2(irho2,j,k,l)
        vsquare = dot_product(prim2(ivel2_x:ivel2_z,j,k,l), prim2(ivel2_x:ivel2_z,j,k,l))
		    epsilon2(j,k,l) = cons2(itau2,j,k,l)/cons2(irho2,j,k,l) - 0.5D0 * vsquare
        factor = MAX(SIGN(1.0D0, epsilon2(j,k,l)), 0.0D0)
        epsilon2(j,k,l) = factor*epsilon2(j,k,l) + (1.0d0 - factor)*eps2_a
	    END DO
    END DO
  END DO
  !$OMP END DO
END IF
!$OMP END PARALLEL

! boundary condition !
call boundary1d_NM (epsilon2,part,even,even,even,even,even,even)
call BOUNDARYP_NM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE