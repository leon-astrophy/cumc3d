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
  DO l = 1, nz
    DO k = 1, ny
      DO j = 1, nx
        bcell(ibx,j,k,l) = 0.5D0*(prim(ibx,j,k,l) + prim(ibx,j-1,k,l))
        bcell(iby,j,k,l) = 0.5D0*(prim(iby,j,k,l) + prim(iby,j,k-1,l))
        bcell(ibz,j,k,l) = 0.5D0*(prim(ibz,j,k,l) + prim(ibz,j,k,l-1))
      END DO
    END DO
  END DO
  !$ACC END PARALLEL
  !$OMP END DO 
ELSEIF(coordinate_flag == 1) THEN
  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(present)
  DO l = 1, nz
    DO k = 1, ny
      DO j = 1, nx
        bcell(ibx,j,k,l) = (xF(j)*prim(ibx,j,k,l) + xF(j-1)*prim(ibx,j-1,k,l))/(xF(j) + xF(j-1))
        bcell(iby,j,k,l) = 0.5D0*(prim(iby,j,k,l) + prim(iby,j,k-1,l))
        bcell(ibz,j,k,l) = 0.5D0*(prim(ibz,j,k,l) + prim(ibz,j,k,l-1))
      END DO
    END DO
  END DO
  !$ACC END PARALLEL
  !$OMP END DO
ELSEIF(coordinate_flag == 2) THEN
  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(present)
  DO l = 1, nz
    DO k = 1, ny
      DO j = 1, nx
        bcell(ibx,j,k,l) = 1.5D0*dx(j)/dx_cb(j)*(xF(j)**2*prim(ibx,j,k,l) + xF(j-1)**2*prim(ibx,j-1,k,l))
        bcell(iby,j,k,l) = 0.5D0*dy(k)/dcose(k)*(sinf(k)*prim(iby,j,k,l) + sinf(k-1)*prim(iby,j,k-1,l))
        bcell(ibz,j,k,l) = 0.5D0*(prim(ibz,j,k,l) + prim(ibz,j,k,l-1))
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
DO l = 0, nz
  DO k = 0, ny
    DO j = 0, nx
      vsquare = dot_product(prim(ivx:ivz,j,k,l), prim(ivx:ivz,j,k,l))
      bsquare = dot_product(bcell(ibx:ibz,j,k,l), bcell(ibx:ibz,j,k,l))
	    cons(irho,j,k,l) = prim(irho,j,k,l)
      cons(ivx:itau-1,j,k,l) = prim(ivx:itau-1,j,k,l)*prim(irho,j,k,l)
	    cons(itau,j,k,l) = prim(irho,j,k,l)*(epsilon(j,k,l) + 0.5D0*vsquare) + 0.5D0*bsquare
      cons(itau+1:ibx-1,j,k,l) = prim(itau+1:ibx-1,j,k,l)*prim(irho,j,k,l)
		  cons(ibx:ibz,j,k,l) = prim(ibx:ibz,j,k,l)
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
! For NM sectors !

!$OMP PARALLEL PRIVATE(vsquare, bsquare)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get cell centered magnetic field !
IF(coordinate_flag == 0) THEN
  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(present)
  DO l = 1, nz
    DO k = 1, ny
      DO j = 1, nx
        bcell(ibx,j,k,l) = 0.5D0*(cons(ibx,j,k,l) + cons(ibx,j-1,k,l))
        bcell(iby,j,k,l) = 0.5D0*(cons(iby,j,k,l) + cons(iby,j,k-1,l))
        bcell(ibz,j,k,l) = 0.5D0*(cons(ibz,j,k,l) + cons(ibz,j,k,l-1))
      END DO
    END DO 
  END DO
  !$ACC END PARALLEL
  !$OMP END DO
ELSEIF(coordinate_flag == 1) THEN
  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(present)
  DO l = 1, nz
    DO k = 1, ny
      DO j = 1, nx
        bcell(ibx,j,k,l) = (xF(j)*cons(ibx,j,k,l) + xF(j-1)*cons(ibx,j-1,k,l))/(xF(j) + xF(j-1))
        bcell(iby,j,k,l) = 0.5D0*(cons(iby,j,k,l) + cons(iby,j,k-1,l))
        bcell(ibz,j,k,l) = 0.5D0*(cons(ibz,j,k,l) + cons(ibz,j,k,l-1))
      END DO
    END DO
  END DO
  !$ACC END PARALLEL
  !$OMP END DO
ELSEIF(coordinate_flag == 2) THEN
  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(present)
  DO l = 1, nz
    DO k = 1, ny
      DO j = 1, nx
        bcell(ibx,j,k,l) = 1.5D0*dx(j)/dx_cb(j)*(xF(j)**2*cons(ibx,j,k,l) + xF(j-1)**2*cons(ibx,j-1,k,l))
        bcell(iby,j,k,l) = 0.5D0*dy(k)/dcose(k)*(sinf(k)*cons(iby,j,k,l) + sinf(k-1)*cons(iby,j,k-1,l))
        bcell(ibz,j,k,l) = 0.5D0*(cons(ibz,j,k,l) + cons(ibz,j,k,l-1))
      END DO
    END DO
  END DO
  !$ACC END PARALLEL
  !$OMP END DO
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do the rest conversation !
!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) PRIVATE(vsquare, bsquare) DEFAULT(present)
DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx
      prim(irho,j,k,l) = cons(irho,j,k,l)
		  prim(ivx:itau-1,j,k,l) = cons(ivx:itau-1,j,k,l)/cons(irho,j,k,l)
      prim(itau+1:ibx-1,j,k,l) = cons(itau+1:ibx-1,j,k,l)/cons(irho,j,k,l)
      prim(ibx:ibz,j,k,l) = cons(ibx:ibz,j,k,l)
      bsquare = dot_product(bcell(ibx:ibz,j,k,l), bcell(ibx:ibz,j,k,l))
      vsquare = dot_product(prim(ivx:ivz,j,k,l), prim(ivx:ivz,j,k,l))
      epsilon(j,k,l) = (cons(itau,j,k,l) - 0.5D0*bsquare)/cons(irho,j,k,l) - 0.5D0 * vsquare
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
