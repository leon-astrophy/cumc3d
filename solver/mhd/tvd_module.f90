!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This module contains subroutines and functions for reconstructing
! interface primitive variables by the 2nd order TVD method
! c.f. Mignone 2014 for more detials 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE TVD_MODULE        
USE DEFINITION
IMPLICIT NONE

! Limiters !
INTEGER, PARAMETER :: minmod = 0
INTEGER, PARAMETER :: vanleer = 1
INTEGER, PARAMETER :: mcentral = 2

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Setup reconstruction weight for TVD
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE TVD_WEIGHT
USE DEFINITION  
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k, l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set up for x, y, and z directional reconstruction !

ALLOCATE(wx(-2:nx+3,1:2))
ALLOCATE(wy(-2:ny+3,1:2))
ALLOCATE(wz(-2:nz+3,1:2))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Assign !

DO j = -1, nx + 2
	wx(j,1) = (dx(j) + dx(j+1))/dx(j)
	wx(j,2) = (dx(j) + dx(j-1))/dx(j)
END DO
DO k = -1, ny + 2
	wy(k,1) = (dy(k) + dy(k+1))/dy(k)
	wy(k,2) = (dy(k) + dy(k-1))/dy(k)
END DO
DO l = -1, nz + 2
	wz(l,1) = (dz(l) + dz(l+1))/dz(l)
	wz(l,2) = (dz(l) + dz(l-1))/dz(l)
END DO

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Do reconstruction using the TVD method along the x-direction, get states at x-interfaces 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE TVD_reconx(limiter)
!$ACC ROUTINE (EOSEPSILON_NM) SEQ
!$ACC ROUTINE (EOSSOUNDSPEED) SEQ
!$ACC ROUTINE (TVD) SEQ
USE RIEMANN_MODULE
USE MHD_MODULE
USE DEFINITION  
IMPLICIT NONE

! Limiter !
INTEGER, INTENT(IN) :: limiter

! Integer !
INTEGER :: j, k, l, i, p

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

! We first interpolate primitive variables for NM !
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE (STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT (PRESENT)
DO l = 1, nz
	DO k = 1, ny
		DO j = 0, nx + 1 

			! NM hydro variables !
			DO i = imin, ibx - 1
				CALL TVD (limiter, wx(j,1:2), prim(i,j-1,k,l), prim(i,j,k,l), prim(i,j+1,k,l), primR(i,j-1,k,l), primL(i,j,k,l))
			END DO

			! Internal energy !
			CALL EOSEPSILON_NM (primR(irho,j-1,k,l), primR(itau,j-1,k,l), epsR(j-1,k,l))
			CALL EOSEPSILON_NM (primL(irho,j,k,l), primL(itau,j,k,l), epsL(j,k,l))
			CALL EOSSOUNDSPEED (primR(itau,j-1,k,l), primR(irho,j-1,k,l), epsR(j-1,k,l), csR(j-1,k,l))
			CALL EOSSOUNDSPEED (primL(itau,j,k,l), primL(irho,j,k,l), epsL(j,k,l), csL(j,k,l))

			! magnetic field !
			CALL TVD (limiter, wx(j,1:2), bcell(iby,j-1,k,l), bcell(iby,j,k,l), bcell(iby,j+1,k,l), primR(iby,j-1,k,l), primL(iby,j,k,l))
			CALL TVD (limiter, wx(j,1:2), bcell(ibz,j-1,k,l), bcell(ibz,j,k,l), bcell(ibz,j+1,k,l), primR(ibz,j-1,k,l), primL(ibz,j,k,l))

			! Special treatment for normal field 
			primR(ibx,j-1,k,l) = prim(ibx,j-1,k,l)
			primL(ibx,j,k,l) = prim(ibx,j,k,l)

		END DO
	END DO
END DO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'tvdx = ', REAL(time_end - time_start) / rate
#endif

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Do reconstruction using the TVD method along the x-direction, get states at x-interfaces 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE TVD_recony(limiter)
!$ACC ROUTINE (EOSEPSILON_NM) SEQ
!$ACC ROUTINE (EOSSOUNDSPEED) SEQ
!$ACC ROUTINE (TVD) SEQ
USE RIEMANN_MODULE
USE MHD_MODULE
USE DEFINITION  
IMPLICIT NONE

! Limiter !
INTEGER, INTENT(IN) :: limiter

! Integer !
INTEGER :: j, k, l, i, p

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

! We first interpolate primitive variables for NM !
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE (STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT (PRESENT)
DO l = 1, nz
	DO k = 0, ny + 1 
		DO j = 1, nx

			! NM hydro variables !
			DO i = imin, ibx - 1
				CALL TVD (limiter, wy(k,1:2), prim(i,j,k-1,l), prim(i,j,k,l), prim(i,j,k+1,l), primR(i,j,k-1,l), primL(i,j,k,l))
			END DO

			! Internal energy !
			CALL EOSEPSILON_NM (primR(irho,j,k-1,l), primR(itau,j,k-1,l), epsR(j,k-1,l))
			CALL EOSEPSILON_NM (primL(irho,j,k,l), primL(itau,j,k,l), epsL(j,k,l))
			CALL EOSSOUNDSPEED (primR(itau,j,k-1,l), primR(irho,j,k-1,l), epsR(j,k-1,l), csR(j,k-1,l))
			CALL EOSSOUNDSPEED (primL(itau,j,k,l), primL(irho,j,k,l), epsL(j,k,l), csL(j,k,l))

			! Magnetic fields !
			CALL TVD (limiter, wy(k,1:2), bcell(ibx,j,k-1,l), bcell(ibx,j,k,l), bcell(ibx,j,k+1,l), primR(ibx,j,k-1,l), primL(ibx,j,k,l))
			CALL TVD (limiter, wy(k,1:2), bcell(ibz,j,k-1,l), bcell(ibz,j,k,l), bcell(ibz,j,k+1,l), primR(ibz,j,k-1,l), primL(ibz,j,k,l))

			! special treatment for formal fields !
			primR(iby,j,k-1,l) = prim(iby,j,k-1,l)
			primL(iby,j,k,l) = prim(iby,j,k,l)	
										 
		END DO
	END DO
END DO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'tvdy = ', REAL(time_end - time_start) / rate
#endif

END SUBROUTINE
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Do reconstruction using the TVD method along the x-direction, get states at x-interfaces 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE TVD_reconz(limiter)
!$ACC ROUTINE (EOSEPSILON_NM) SEQ
!$ACC ROUTINE (EOSSOUNDSPEED) SEQ
!$ACC ROUTINE (TVD) SEQ
USE RIEMANN_MODULE
USE MHD_MODULE
USE DEFINITION  
IMPLICIT NONE

! Limiter !
INTEGER, INTENT(IN) :: limiter

! Integer !
INTEGER :: j, k, l, i, p

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

! We first interpolate primitive variables for NM !
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE (STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT (PRESENT)
DO l = 0, nz + 1 
	DO k = 1, ny
		DO j = 1, nx

			! NM hydro variables !
			DO i = imin, ibx - 1
				CALL TVD (limiter, wz(l,1:2), prim(i,j,k,l-1), prim(i,j,k,l), prim(i,j,k,l+1), primR(i,j,k,l-1), primL(i,j,k,l))
			END DO

			! Internal energy !
			CALL EOSEPSILON_NM (primR(irho,j,k,l-1), primR(itau,j,k,l-1), epsR(j,k,l-1))
			CALL EOSEPSILON_NM (primL(irho,j,k,l), primL(itau,j,k,l), epsL(j,k,l))
			CALL EOSSOUNDSPEED (primR(itau,j,k,l-1), primR(irho,j,k,l-1), epsR(j,k,l-1), csR(j,k,l-1))
			CALL EOSSOUNDSPEED (primL(itau,j,k,l), primL(irho,j,k,l), epsL(j,k,l), csL(j,k,l))
	
			! Magnetic field !
			CALL TVD (limiter, wz(l,1:2), bcell(ibx,j,k,l-1), bcell(ibx,j,k,l), bcell(ibx,j,k,l+1), primR(ibx,j,k,l-1), primL(ibx,j,k,l))
			CALL TVD (limiter, wz(l,1:2), bcell(iby,j,k,l-1), bcell(iby,j,k,l), bcell(iby,j,k,l+1), primR(iby,j,k,l-1), primL(iby,j,k,l))

			! special treatment for normal field !
			primR(ibz,j,k,l-1) = prim(ibz,j,k,l-1)
			primL(ibz,j,k,l) = prim(ibz,j,k,l)										

		END DO
	END DO
END DO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'tvdz = ', REAL(time_end - time_start) / rate
#endif

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The modified van-leer (VL) slope limiter for getting interface primitive variables 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE TVD (lim_in, wgt, vm, vc, vp, vm_out, vp_out)
!$ACC ROUTINE SEQ 
USE DEFINITION
IMPLICIT NONE

! Input conservative variables !
INTEGER, INTENT(IN) :: lim_in
REAL*8, INTENT (IN) :: vm, vc, vp

! The reconstructed states at cell boundary !
REAL*8, INTENT (OUT) :: vm_out, vp_out

! Grid size !
REAL*8, DIMENSION(1:2) :: wgt

! Integer !
INTEGER :: i, j, k

! Real !
REAL*8 :: dm1, dc, dp1
REAL*8 :: cf, cb, dq2
REAL*8 :: dqb, dqf
REAL*8 :: slope, v_tvd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign !
cf = wgt(1)
cb = wgt(2)

! assign !
dqf = 2.0d0*(vp - vc)/cf
dqb = 2.0d0*(vc - vm)/cb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Choose by limiter 

! Assign slope !
IF(lim_in == minmod) THEN
	slope = 0.5d0*(SIGN(1.0d0, dqf) + SIGN(1.0d0, dqb))*MIN(ABS(dqf), ABS(dqb))
ELSEIF(lim_in == vanleer) THEN
	dq2 = dqf*dqb
	IF (dq2 <= 0.0D0) THEN
		slope = 0.0
	ELSE
		slope = (dq2*(cf*dqb + cb*dqf)/(dqb*dqb + dqf*dqf + dq2*(cf + cb - 2.0)))
	END IF
ELSEIF(lim_in == mcentral) THEN
	IF(dqf == 0.0D0) THEN
		slope = 0.0D0
	ELSE
		v_tvd = dqb/dqf
		slope = MAX(0.0D0, MIN(0.5D0*(1.0D0 + v_tvd), cf, cb*v_tvd))
		slope = dqf*slope
	END IF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Output !
vp_out = vc + 0.5d0*slope
vm_out = vc - 0.5d0*slope

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE
