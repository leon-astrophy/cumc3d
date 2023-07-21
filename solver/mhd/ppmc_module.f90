!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the PPM (Piecewise Parabolic Method) Module that reconstruct   
! interface values of primitive variables. Three different appoarch 
! can be chosen, and they differ by how to handle monotone conditions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE PPMC_MODULE
USE DEFINITION
IMPLICIT NONE

contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do reconstruction using PPM interpolation, with the original Colella algorithm
! Reconstruct along the x-direction, and get states at the x-interfaces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PPMC_reconx
!$ACC ROUTINE (EOSEPSILON_NM) SEQ
!$ACC ROUTINE (EOSSOUNDSPEED) SEQ
!$ACC ROUTINE (PPMC) SEQ
USE RIEMANN_MODULE
USE MHD_MODULE
USE DEFINITION  
IMPLICIT NONE

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
DO l = nz_min_2 - 1, nz_part_2 + 1
	DO k = ny_min_2 - 1, ny_part_2 + 1
		DO j = nx_min_2 - 1, nx_part_2 + 1

			! NM hydro variables !
			DO i = imin2, ibx - 1
				CALL PPMC (dx2(j-2), dx2(j-1), dx2(j), dx2(j+1), dx2(j+2), & 
									 prim2(i,j-2,k,l), prim2(i,j-1,k,l), prim2(i,j,k,l), prim2(i,j+1,k,l), prim2(i,j+2,k,l), & 
									 primR2(i,j-1,k,l), primL2(i,j,k,l))
			END DO

			! Internal energy !
			CALL EOSEPSILON_NM (primR2(irho2,j-1,k,l), primR2(itau2,j-1,k,l), eps2R(j-1,k,l))
			CALL EOSEPSILON_NM (primL2(irho2,j,k,l), primL2(itau2,j,k,l), eps2L(j,k,l))
			CALL EOSSOUNDSPEED (primR2(itau2,j-1,k,l), primR2(irho2,j-1,k,l), eps2R(j-1,k,l), cs2R(j-1,k,l))
			CALL EOSSOUNDSPEED (primL2(itau2,j,k,l), primL2(irho2,j,k,l), eps2L(j,k,l), cs2L(j,k,l))

			! magnetic field !
			CALL PPMC (dx2(j-2), dx2(j-1), dx2(j), dx2(j+1), dx2(j+2), & 
								 bcell(iby,j-2,k,l), bcell(iby,j-1,k,l), bcell(iby,j,k,l), bcell(iby,j+1,k,l), bcell(iby,j+2,k,l), & 
								 primR2(iby,j-1,k,l), primL2(iby,j,k,l))
			CALL PPMC (dx2(j-2), dx2(j-1), dx2(j), dx2(j+1), dx2(j+2), & 
								 bcell(ibz,j-2,k,l), bcell(ibz,j-1,k,l), bcell(ibz,j,k,l), bcell(ibz,j+1,k,l), bcell(ibz,j+2,k,l), & 
								 primR2(ibz,j-1,k,l), primL2(ibz,j,k,l))

			! Special treatment for normal field 
			primR2(ibx,j-1,k,l) = prim2(ibx,j-1,k,l)
			primL2(ibx,j,k,l) = prim2(ibx,j,k,l)

		END DO
	END DO
END DO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'ppmcx = ', REAL(time_end - time_start) / rate
#endif

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do reconstruction using PPM interpolation, with the original Colella algorithm
! Reconstruct along the y-direction, and get states at the y-interfaces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PPMC_recony
!$ACC ROUTINE (EOSEPSILON_NM) SEQ
!$ACC ROUTINE (EOSSOUNDSPEED) SEQ
!$ACC ROUTINE (PPMC) SEQ
USE RIEMANN_MODULE
USE MHD_MODULE
USE DEFINITION  
IMPLICIT NONE

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
DO l = nz_min_2 - 1, nz_part_2 + 1
	DO k = ny_min_2 - 1, ny_part_2 + 1
		DO j = nx_min_2 - 1, nx_part_2 + 1

			! NM hydro variables !
			DO i = imin2, ibx - 1
				CALL PPMC (dy2(k-2), dy2(k-1), dy2(k), dy2(k+1), dy2(k+2), & 
									 prim2(i,j,k-2,l), prim2(i,j,k-1,l), prim2(i,j,k,l), prim2(i,j,k+1,l), prim2(i,j,k+2,l), & 
									 primR2(i,j,k-1,l), primL2(i,j,k,l))
			END DO

			! Internal energy !
			CALL EOSEPSILON_NM (primR2(irho2,j,k-1,l), primR2(itau2,j,k-1,l), eps2R(j,k-1,l))
			CALL EOSEPSILON_NM (primL2(irho2,j,k,l), primL2(itau2,j,k,l), eps2L(j,k,l))
			CALL EOSSOUNDSPEED (primR2(itau2,j,k-1,l), primR2(irho2,j,k-1,l), eps2R(j,k-1,l), cs2R(j,k-1,l))
			CALL EOSSOUNDSPEED (primL2(itau2,j,k,l), primL2(irho2,j,k,l), eps2L(j,k,l), cs2L(j,k,l))

			! Magnetic fields !
			CALL PPMC (dy2(k-2), dy2(k-1), dy2(k), dy2(k+1), dy2(k+2), & 
								 bcell(ibx,j,k-2,l), bcell(ibx,j,k-1,l), bcell(ibx,j,k,l), bcell(ibx,j,k+1,l), bcell(ibx,j,k+2,l), & 
								 primR2(ibx,j,k-1,l), primL2(ibx,j,k,l))
			CALL PPMC (dy2(k-2), dy2(k-1), dy2(k), dy2(k+1), dy2(k+2), & 
								 bcell(ibz,j,k-2,l), bcell(ibz,j,k-1,l), bcell(ibz,j,k,l), bcell(ibz,j,k+1,l), bcell(ibz,j,k+2,l), & 
								 primR2(ibz,j,k-1,l), primL2(ibz,j,k,l))

			! special treatment for formal fields !
			primR2(iby,j,k-1,l) = prim2(iby,j,k-1,l)
			primL2(iby,j,k,l) = prim2(iby,j,k,l)	
										 
		END DO
	END DO
END DO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'ppmcy = ', REAL(time_end - time_start) / rate
#endif

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do reconstruction using PPM interpolation, with the original Colella algorithm
! Reconstruct along the z-direction, and get states at the z-interfaces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PPMC_reconz
!$ACC ROUTINE (EOSEPSILON_NM) SEQ
!$ACC ROUTINE (EOSSOUNDSPEED) SEQ
!$ACC ROUTINE (PPMC) SEQ
USE RIEMANN_MODULE
USE MHD_MODULE
USE DEFINITION  
IMPLICIT NONE

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
DO l = nz_min_2 - 1, nz_part_2 + 1
	DO k = ny_min_2 - 1, ny_part_2 + 1
		DO j = nx_min_2 - 1, nx_part_2 + 1

			! NM hydro variables !
			DO i = imin2, ibx - 1
				CALL PPMC (dz2(l-2), dz2(l-1), dz2(l), dz2(l+1), dz2(l+2), & 
									 prim2(i,j,k,l-2), prim2(i,j,k,l-1), prim2(i,j,k,l), prim2(i,j,k,l+1), prim2(i,j,k,l+2), & 
									 primR2(i,j,k,l-1), primL2(i,j,k,l))
			END DO

			! Internal energy !
			CALL EOSEPSILON_NM (primR2(irho2,j,k,l-1), primR2(itau2,j,k,l-1), eps2R(j,k,l-1))
			CALL EOSEPSILON_NM (primL2(irho2,j,k,l), primL2(itau2,j,k,l), eps2L(j,k,l))
			CALL EOSSOUNDSPEED (primR2(itau2,j,k,l-1), primR2(irho2,j,k,l-1), eps2R(j,k,l-1), cs2R(j,k,l-1))
			CALL EOSSOUNDSPEED (primL2(itau2,j,k,l), primL2(irho2,j,k,l), eps2L(j,k,l), cs2L(j,k,l))
	
			! Magnetic field !
			CALL PPMC (dz2(l-2), dz2(l-1), dz2(l), dz2(l+1), dz2(l+2), & 
								 bcell(ibx,j,k,l-2), bcell(ibx,j,k,l-1), bcell(ibx,j,k,l), bcell(ibx,j,k,l+1), bcell(ibx,j,k,l+2), & 
								 primR2(ibx,j,k,l-1), primL2(ibx,j,k,l))
			CALL PPMC (dz2(l-2), dz2(l-1), dz2(l), dz2(l+1), dz2(l+2), & 
								 bcell(iby,j,k,l-2), bcell(iby,j,k,l-1), bcell(iby,j,k,l), bcell(iby,j,k,l+1), bcell(iby,j,k,l+2), & 
								 primR2(iby,j,k,l-1), primL2(iby,j,k,l))

			! special treatment for normal field !
			primR2(ibz,j,k,l-1) = prim2(ibz,j,k,l-1)
			primL2(ibz,j,k,l) = prim2(ibz,j,k,l)										

		END DO
	END DO
END DO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'ppmcz = ', REAL(time_end - time_start) / rate
#endif

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate cell average values to interface values using PPM 
! Assumed non-uniform gridding. Using the original PPM algorithm
! See Colella 1984. No steepening or flattening is performed.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PPMC (dm2, dm1, dc, dp1, dp2, vm2, vm1, vc, vp1, vp2, vm_out, vp_out)
!$ACC ROUTINE SEQ 
USE DEFINITION
IMPLICIT NONE

! The input into the subroutine, including conservative variable and input flux function !
REAL*8, INTENT (IN) :: dm2, dm1, dc, dp1, dp2
REAL*8, INTENT (IN) :: vm2, vm1, vc, vp1, vp2

! The output of the subroutine, the flux at cell boundary !
REAL*8, INTENT (OUT) :: vp_out, vm_out

! Temporal variables !
REAL*8 :: vl, vr
REAL*8 :: dmp1, dmc, dmm1
REAL*8 :: deltap1, deltac, deltam1
REAL*8 :: a1R, a2R, deltaXR, z1R, z2R
REAL*8 :: a1L, a2L, deltaXL, z1L, z2L
REAL*8 :: cp1, cp2, cp3
REAL*8 :: cc1, cc2, cc3
REAL*8 :: cm1, cm2, cm3
REAL*8 :: condition

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Get coefficient for interpolants!
a1R = dc/(dc + dp1)
a1L = dm1/(dm1 + dc)
a2R = 2.0D0*dp1*dc/(dp1 + dc)
a2L = 2.0D0*dc*dm1/(dc + dm1)
deltaXR = dm1 + dc + dp1 + dp2
deltaXL = dm2 + dm1 + dc + dp1
z1R = (dm1 + dc)/(2.0D0*dc + dp1)
z1L = (dm2 + dm1)/(2.0D0*dm1 + dc)
z2R = (dp2 + dp1)/(2.0D0*dp1 + dc)
z2L = (dp1 + dc)/(2.0D0*dc + dm1)

! For slope estimations !
cc1 = dc/(dm1 + dc + dp1)
cp1 = dp1/(dc + dp1 + dp2)
cm1 = dm1/(dm2 + dm1 + dc)
cc2 = (2.0D0*dm1 + dc)/(dp1 + dc)
cp2 = (2.0D0*dc + dp1)/(dp2 + dp1)
cm2 = (2.0D0*dm2 + dm1)/(dc + dm1)
cc3 = (dc + 2.0D0*dp1)/(dm1 + dc)
cp3 = (dp1 + 2.0D0*dp2)/(dc + dp1)
cm3 = (dm1 + 2.0D0*dc)/(dm2 + dm1)

! Get slopes !
deltap1 = cp1*(cp2*(vp2 - vp1) + cp3*(vp1 - vc))
deltac = cc1*(cc2*(vp1 - vc) + cc3*(vc - vm1))
deltam1 = cm1*(cm2*(vc - vm1) + cm3*(vm1 - vm2))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Van Leer limiter !

condition = (vp1 - vc)*(vc - vm1)
IF(condition > 0.0D0) THEN
	dmc = min(abs(deltac), 2.0D0*abs(vc - vm1), 2.0D0*abs(vc - vp1))*sign(1.0D0,deltac)
ELSE
	dmc = 0.0D0
END IF

condition = (vp2 - vp1)*(vp1 - vc)
IF(condition > 0.0D0) THEN
	dmp1 = min(abs(deltap1), 2.0D0*abs(vp1 - vc), 2.0D0*abs(vp1 - vp2))*sign(1.0D0,deltap1)
ELSE
	dmp1 = 0.0D0
END IF

condition = (vc - vm1)*(vm1 - vm2)
IF(condition > 0.0D0) THEN
	dmm1 = min(abs(deltam1), 2.0D0*abs(vm1 - vm2), 2.0D0*abs(vm1 - vc))*sign(1.0D0,deltam1)
ELSE
	dmm1 = 0.0D0
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Get the interpolant !
vp_out = vc + a1R*(vp1 - vc) + (a2R*(z1R - z2R)*(vp1 - vc) - dc*z1R*dmp1 + dp1*z2R*dmc)/deltaXR
vm_out = vm1 + a1L*(vc - vm1) + (a2L*(z1L - z2L)*(vc - vm1) - dm1*z1L*dmc + dc*z2L*dmm1)/deltaXL

! backup !
vl = vm_out
vr = vp_out

! Check for extremum conditions !
IF((vr - vc)*(vc - vl) <= 0.0D0) THEN
	vp_out = vc
	vm_out = vc
ELSE
	condition = (vr - vl)*(vl - 3.0D0*vc + 2.0D0*vr)
	IF(condition < 0.0D0) THEN
		vm_out = 3.0d0*vc - 2.0D0*vr
	END IF
	condition = (vp_out - vm_out)*(3.0D0*vc - 2.0d0*vl - vr)
	IF(condition < 0.0D0) THEN
		vp_out = 3.0d0*vc - 2.0D0*vl
	END IF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

END MODULE
