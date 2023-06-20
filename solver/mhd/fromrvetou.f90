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

! Check timing with or without openmp
INTEGER :: time_start, time_end
INTEGER :: cr
REAL*8 :: rate
CALL system_clock(count_rate=cr)
rate = REAL(cr)

CALL system_clock(time_start)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NM Sector !

! Convert NM hydro
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(vsquare, bsquare)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) PRIVATE(vsquare, bsquare) DEFAULT(present)
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
!$ACC END PARALLEL
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL system_clock(time_end)
#ifdef DEBUG
WRITE(*,*) 'p2u = ', REAL(time_end - time_start) / rate
#endif

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
REAL*8 :: diff, factor

! Check timing with or without openmp
INTEGER :: time_start, time_end
INTEGER :: cr
REAL*8 :: rate
CALL system_clock(count_rate=cr)
rate = REAL(cr)

CALL system_clock(time_start)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For NM sectors !

! Convert the NM hydro
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(vsquare, bsquare, diff, factor, et_2, bige_2)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) PRIVATE(vsquare, bsquare, diff, factor, et_2, bige_2) DEFAULT(present)
DO l = nz_min_2, nz_part_2
  DO k = ny_min_2, ny_part_2
    DO j = nx_min_2, nx_part_2
		  prim2(imin2:imax2,j,k,l) = cons2(imin2:imax2,j,k,l)/cons2(irho2,j,k,l)
	    prim2(irho2,j,k,l) = cons2(irho2,j,k,l)
      prim2(ibx:ibz,j,k,l) = cons2(ibx:ibz,j,k,l)
      bsquare = dot_product(prim2(ibx:ibz,j,k,l), prim2(ibx:ibz,j,k,l))
      cons2(itau2,j,k,l) = cons2(itau2,j,k,l) - 0.5D0*bsquare
      vsquare = dot_product(prim2(ivel2_x:ivel2_z,j,k,l), prim2(ivel2_x:ivel2_z,j,k,l))
      epsilon2(j,k,l) = cons2(itau2,j,k,l)/cons2(irho2,j,k,l) - 0.5D0 * vsquare
      factor = MAX(SIGN(1.0D0, epsilon2(j,k,l)), 0.0D0)
      epsilon2(j,k,l) = factor*epsilon2(j,k,l) + (1.0d0 - factor)*eps2_a
    END DO
  END DO
END DO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! boundary condition !
call boundary1d_NM (epsilon2, part, even, even, even, even, even, even)
call BOUNDARYP_NM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL system_clock(time_end)
#ifdef DEBUG
WRITE(*,*) 'u2p = ', REAL(time_end - time_start) / rate
#endif

END SUBROUTINE