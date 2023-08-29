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
#ifdef DEBUG
INTEGER :: time_start, time_end
INTEGER :: cr
REAL*8 :: rate
CALL system_clock(count_rate=cr)
rate = REAL(cr)
CALL system_clock(time_start)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NM Sector !

!$OMP PARALLEL PRIVATE(vsquare, bsquare)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get cell centered magnetic field !
IF(coordinate_flag == 0) THEN
  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(present)
  DO l = nz_min_2, nz_part_2
    DO k = ny_min_2, ny_part_2
      DO j = nx_min_2, nx_part_2
        bcell(ibx,j,k,l) = 0.5D0*(prim2(ibx,j,k,l) + prim2(ibx,j-1,k,l))
        bcell(iby,j,k,l) = 0.5D0*(prim2(iby,j,k,l) + prim2(iby,j,k-1,l))
        bcell(ibz,j,k,l) = 0.5D0*(prim2(ibz,j,k,l) + prim2(ibz,j,k,l-1))
      END DO
    END DO
  END DO
  !$ACC END PARALLEL
  !$OMP END DO
ELSEIF(coordinate_flag == 1) THEN
  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(present)
  DO l = nz_min_2, nz_part_2
    DO k = ny_min_2, ny_part_2
      DO j = nx_min_2, nx_part_2
        bcell(ibx,j,k,l) = ((xF2(j) - x2cen(j))*prim2(ibx,j,k,l) + (x2cen(j) - xF2(j-1))*prim2(ibx,j-1,k,l))/(dx2(j))
        !(xF2(j)*prim2(ibx,j,k,l) + xF2(j-1)*prim2(ibx,j-1,k,l))/(xF2(j) + xF2(j-1))
        bcell(iby,j,k,l) = ((yF2(k) - y2cen(k))*prim2(iby,j,k,l) + (y2cen(k) - yF2(k-1))*prim2(iby,j,k-1,l))/(dy2(k))
        !0.5D0*(prim2(iby,j,k,l) + prim2(iby,j,k-1,l))
        bcell(ibz,j,k,l) = 0.5D0*(prim2(ibz,j,k,l) + prim2(ibz,j,k,l-1))
      END DO
    END DO
  END DO
  !$ACC END PARALLEL
  !$OMP END DO
ELSEIF(coordinate_flag == 2) THEN
  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(present)
  DO l = nz_min_2, nz_part_2
    DO k = ny_min_2, ny_part_2
      DO j = nx_min_2, nx_part_2
        bcell(ibx,j,k,l) = ((xF2(j) - x2cen(j))*prim2(ibx,j,k,l) + (x2cen(j) - xF2(j-1))*prim2(ibx,j-1,k,l))/(dx2(j))
        !1.5D0*dx2(j)/dx2_cb(j)*(xF2(j)**2*prim2(ibx,j,k,l) + xF2(j-1)**2*prim2(ibx,j-1,k,l))
        bcell(iby,j,k,l) = ((yF2(k) - y2cen(k))*prim2(iby,j,k,l) + (y2cen(k) - yF2(k-1))*prim2(iby,j,k-1,l))/(dy2(k))
        !0.5D0*dy2(k)/dcos2(k)*(sin2f(k)*prim2(iby,j,k,l) + sin2f(k-1)*prim2(iby,j,k-1,l))
        bcell(ibz,j,k,l) = 0.5D0*(prim2(ibz,j,k,l) + prim2(ibz,j,k,l-1))
      END DO
    END DO
  END DO
  !$ACC END PARALLEL
  !$OMP END DO
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Convert NM hydro
!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(present) PRIVATE(vsquare, bsquare)
DO l = nz_min_2 - 1, nz_part_2
  DO k = ny_min_2 - 1, ny_part_2
    DO j = nx_min_2 - 1, nx_part_2
      vsquare = dot_product(prim2(ivel2_x:ivel2_z,j,k,l), prim2(ivel2_x:ivel2_z,j,k,l))
      bsquare = dot_product(bcell(ibx:ibz,j,k,l), bcell(ibx:ibz,j,k,l))
		  cons2(imin2:ibx-1,j,k,l) = prim2(imin2:ibx-1,j,k,l)*prim2(irho2,j,k,l)
	    cons2(irho2,j,k,l) = prim2(irho2,j,k,l)
	    cons2(itau2,j,k,l) = prim2(irho2,j,k,l)*(epsilon2(j,k,l) + 0.5D0*vsquare) + 0.5D0*bsquare
		  cons2(ibx:ibz,j,k,l) = prim2(ibx:ibz,j,k,l)
	  END DO
  END DO
END DO
!$ACC END PARALLEL
!$OMP END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$OMP END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'p2u = ', REAL(time_end - time_start) / rate
#endif

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Convert back to primitive variables !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FROMUTORVE
USE DEFINITION
USE MHD_MODULE
use ieee_arithmetic
IMPLICIT NONE

! Dummy variables
INTEGER :: i, j, k, l

! Squared velocity !
REAL*8 :: vsquare, bsquare

! For epsilon !
REAL*8 :: factor

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
! For NM sectors !

!$OMP PARALLEL PRIVATE(vsquare, bsquare, factor)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Convert the NM hydro
!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(present)
DO l = nz_min_2 - 1, nz_part_2
  DO k = ny_min_2 - 1, ny_part_2
    DO j = nx_min_2 - 1, nx_part_2
		  prim2(imin2:ibx-1,j,k,l) = cons2(imin2:ibx-1,j,k,l)/cons2(irho2,j,k,l)
      prim2(ibx:ibz,j,k,l) = cons2(ibx:ibz,j,k,l)
    END DO
  END DO
END DO
!$ACC END PARALLEL
!$OMP END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get cell centered magnetic field !
IF(coordinate_flag == 0) THEN
  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(present)
  DO l = nz_min_2, nz_part_2
    DO k = ny_min_2, ny_part_2
      DO j = nx_min_2, nx_part_2
        bcell(ibx,j,k,l) = 0.5D0*(prim2(ibx,j,k,l) + prim2(ibx,j-1,k,l))
        bcell(iby,j,k,l) = 0.5D0*(prim2(iby,j,k,l) + prim2(iby,j,k-1,l))
        bcell(ibz,j,k,l) = 0.5D0*(prim2(ibz,j,k,l) + prim2(ibz,j,k,l-1))
      END DO
    END DO 
  END DO
  !$ACC END PARALLEL
  !$OMP END DO
ELSEIF(coordinate_flag == 1) THEN
  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(present)
  DO l = nz_min_2, nz_part_2
    DO k = ny_min_2, ny_part_2
      DO j = nx_min_2, nx_part_2
        bcell(ibx,j,k,l) = ((xF2(j) - x2cen(j))*prim2(ibx,j,k,l) + (x2cen(j) - xF2(j-1))*prim2(ibx,j-1,k,l))/(dx2(j))
        !(xF2(j)*prim2(ibx,j,k,l) + xF2(j-1)*prim2(ibx,j-1,k,l))/(xF2(j) + xF2(j-1))
        bcell(iby,j,k,l) = ((yF2(k) - y2cen(k))*prim2(iby,j,k,l) + (y2cen(k) - yF2(k-1))*prim2(iby,j,k-1,l))/(dy2(k))
        !0.5D0*(prim2(iby,j,k,l) + prim2(iby,j,k-1,l))
        bcell(ibz,j,k,l) = 0.5D0*(prim2(ibz,j,k,l) + prim2(ibz,j,k,l-1))
      END DO
    END DO
  END DO
  !$ACC END PARALLEL
  !$OMP END DO
ELSEIF(coordinate_flag == 2) THEN
  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(present)
  DO l = nz_min_2, nz_part_2
    DO k = ny_min_2, ny_part_2
      DO j = nx_min_2, nx_part_2
        bcell(ibx,j,k,l) = ((xF2(j) - x2cen(j))*prim2(ibx,j,k,l) + (x2cen(j) - xF2(j-1))*prim2(ibx,j-1,k,l))/(dx2(j))
        !1.5D0*dx2(j)/dx2_cb(j)*(xF2(j)**2*prim2(ibx,j,k,l) + xF2(j-1)**2*prim2(ibx,j-1,k,l))
        bcell(iby,j,k,l) = ((yF2(k) - y2cen(k))*prim2(iby,j,k,l) + (y2cen(k) - yF2(k-1))*prim2(iby,j,k-1,l))/(dy2(k))
        !0.5D0*dy2(k)/dcos2(k)*(sin2f(k)*prim2(iby,j,k,l) + sin2f(k-1)*prim2(iby,j,k-1,l))
        bcell(ibz,j,k,l) = 0.5D0*(prim2(ibz,j,k,l) + prim2(ibz,j,k,l-1))
      END DO
    END DO
  END DO
  !$ACC END PARALLEL
  !$OMP END DO
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do the rest conversation !
!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) PRIVATE(vsquare, bsquare, factor) DEFAULT(present)
DO l = nz_min_2, nz_part_2
  DO k = ny_min_2, ny_part_2
    DO j = nx_min_2, nx_part_2
      prim2(irho2,j,k,l) = cons2(irho2,j,k,l)
      bsquare = dot_product(bcell(ibx:ibz,j,k,l), bcell(ibx:ibz,j,k,l))
      vsquare = dot_product(prim2(ivel2_x:ivel2_z,j,k,l), prim2(ivel2_x:ivel2_z,j,k,l))
      epsilon2(j,k,l) = (cons2(itau2,j,k,l) - 0.5D0*bsquare)/cons2(irho2,j,k,l) - 0.5D0 * vsquare
      factor = MAX(SIGN(1.0D0, epsilon2(j,k,l)), 0.0D0)
      epsilon2(j,k,l) = factor*epsilon2(j,k,l) + (1.0d0 - factor)*eps2_a
    END DO
  END DO
END DO
!$ACC END PARALLEL
!$OMP END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$OMP END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'u2p = ', REAL(time_end - time_start) / rate
#endif

END SUBROUTINE
