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
!
! Set up PPM reconstruction weight 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PPM_WEIGHT
USE DEFINITION  
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k, l

! Reak !
REAL*8 :: dm2, dm1, dc, dp1, dp2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set up for x, y, and z directional reconstruction !

ALLOCATE(wx(-2:nx+3,14))
ALLOCATE(wy(-2:ny+3,14))
ALLOCATE(wz(-2:nz+3,14))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Assign  reconstruction weight for the x-direction !
DO j = 0, nx + 1

	! Grid !
	dm2 = dx(j-2)
	dm1 = dx(j-1)
	dc = dx(j)
	dp1 = dx(j+1)
	dp2 = dx(j+2)

	! Assign weight !
	CALL PPM_INTERPOLANT(dm2, dm1, dc, dp1, dp2, wx(j,1:14))

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Assign  reconstruction weight for the y-direction !
DO j = 0, ny + 1

	! Grid !
	dm2 = dy(j-2)
	dm1 = dy(j-1)
	dc = dy(j)
	dp1 = dy(j+1)
	dp2 = dy(j+2)

	! Assign weight !
	CALL PPM_INTERPOLANT(dm2, dm1, dc, dp1, dp2, wy(j,1:14))

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Assign  reconstruction weight for the z-direction !
DO j = 0, nz + 1

	! Grid !
	dm2 = dz(j-2)
	dm1 = dz(j-1)
	dc = dz(j)
	dp1 = dz(j+1)
	dp2 = dz(j+2)

	! Assign weight !
	CALL PPM_INTERPOLANT(dm2, dm1, dc, dp1, dp2, wz(j,1:14))

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Get PPM interpolant weight from the uneven grid size dx
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PPM_INTERPOLANT(dm2, dm1, dc, dp1, dp2, wgt_out)
USE DEFINITION
IMPLICIT NONE

! Input !
REAL*8, INTENT(IN) :: dm2, dm1, dc, dp1, dp2
REAL*8, DIMENSION(1:14) :: wgt_out

! Real !
REAL*8 :: a1R, aR, deltaXR, z1R, zR
REAL*8 :: a1L, aL, deltaXL, z1L, zL
REAL*8 :: cp1, cp2, cp3
REAL*8 :: cc1, cc2, cc3
REAL*8 :: cm1, cm2, cm3

! Get coefficient for interpolants!
a1R = dc/(dc + dp1)
a1L = dm1/(dm1 + dc)
aR = 2.0D0*dp1*dc/(dp1 + dc)
aL = 2.0D0*dc*dm1/(dc + dm1)
deltaXR = dm1 + dc + dp1 + dp2
deltaXL = dm2 + dm1 + dc + dp1
z1R = (dm1 + dc)/(2.0D0*dc + dp1)
z1L = (dm2 + dm1)/(2.0D0*dm1 + dc)
zR = (dp2 + dp1)/(2.0D0*dp1 + dc)
zL = (dp1 + dc)/(2.0D0*dc + dm1)

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

! Assign weight !
wgt_out(1) = a1L
wgt_out(2) = aL*(z1L - zL)/deltaXL
wgt_out(3) = (-dm1*z1L)/deltaXL
wgt_out(4) = dc*zL/deltaXL
wgt_out(5) = a1R
wgt_out(6) = aR*(z1R - zR)/deltaXR
wgt_out(7) = (-dc*z1R)/deltaXR
wgt_out(8) = dp1*zR/deltaXR
wgt_out(9) = cm1*cm2
wgt_out(10) = cm1*cm3
wgt_out(11) = cc1*cc2
wgt_out(12) = cc1*cc3
wgt_out(13) = cp1*cp2
wgt_out(14) = cp1*cp3

END SUBROUTINE

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
DO l = 1, nz
	DO k = 1, ny
		DO j = 0, nx + 1 

			! NM hydro variables !
			DO i = imin, ibx - 1
				CALL PPMC (wx(j,1:14), prim(i,j-2,k,l), prim(i,j-1,k,l), prim(i,j,k,l), prim(i,j+1,k,l), prim(i,j+2,k,l), primR(i,j-1,k,l), primL(i,j,k,l))
			END DO

			! Internal energy !
			CALL EOSEPSILON_NM (primR(irho,j-1,k,l), primR(itau,j-1,k,l), epsR(j-1,k,l))
			CALL EOSEPSILON_NM (primL(irho,j,k,l), primL(itau,j,k,l), epsL(j,k,l))
			CALL EOSSOUNDSPEED (primR(itau,j-1,k,l), primR(irho,j-1,k,l), epsR(j-1,k,l), csR(j-1,k,l))
			CALL EOSSOUNDSPEED (primL(itau,j,k,l), primL(irho,j,k,l), epsL(j,k,l), csL(j,k,l))

			! magnetic field !
			CALL PPMC (wx(j,1:14), bcell(iby,j-2,k,l), bcell(iby,j-1,k,l), bcell(iby,j,k,l), bcell(iby,j+1,k,l), bcell(iby,j+2,k,l), primR(iby,j-1,k,l), primL(iby,j,k,l))
			CALL PPMC (wx(j,1:14), bcell(ibz,j-2,k,l), bcell(ibz,j-1,k,l), bcell(ibz,j,k,l), bcell(ibz,j+1,k,l), bcell(ibz,j+2,k,l), primR(ibz,j-1,k,l), primL(ibz,j,k,l))

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
DO l = 1, nz
	DO k = 0, ny + 1 
		DO j = 1, nx

			! NM hydro variables !
			DO i = imin, ibx - 1
				CALL PPMC (wy(k,1:14), prim(i,j,k-2,l), prim(i,j,k-1,l), prim(i,j,k,l), prim(i,j,k+1,l), prim(i,j,k+2,l), primR(i,j,k-1,l), primL(i,j,k,l))
			END DO

			! Internal energy !
			CALL EOSEPSILON_NM (primR(irho,j,k-1,l), primR(itau,j,k-1,l), epsR(j,k-1,l))
			CALL EOSEPSILON_NM (primL(irho,j,k,l), primL(itau,j,k,l), epsL(j,k,l))
			CALL EOSSOUNDSPEED (primR(itau,j,k-1,l), primR(irho,j,k-1,l), epsR(j,k-1,l), csR(j,k-1,l))
			CALL EOSSOUNDSPEED (primL(itau,j,k,l), primL(irho,j,k,l), epsL(j,k,l), csL(j,k,l))

			! Magnetic fields !
			CALL PPMC (wy(k,1:14), bcell(ibx,j,k-2,l), bcell(ibx,j,k-1,l), bcell(ibx,j,k,l), bcell(ibx,j,k+1,l), bcell(ibx,j,k+2,l), primR(ibx,j,k-1,l), primL(ibx,j,k,l))
			CALL PPMC (wy(k,1:14), bcell(ibz,j,k-2,l), bcell(ibz,j,k-1,l), bcell(ibz,j,k,l), bcell(ibz,j,k+1,l), bcell(ibz,j,k+2,l), primR(ibz,j,k-1,l), primL(ibz,j,k,l))

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
DO l = 0, nz + 1 
	DO k = 1, ny
		DO j = 1, nx

			! NM hydro variables !
			DO i = imin, ibx - 1
				CALL PPMC (wz(l,1:14), prim(i,j,k,l-2), prim(i,j,k,l-1), prim(i,j,k,l), prim(i,j,k,l+1), prim(i,j,k,l+2), primR(i,j,k,l-1), primL(i,j,k,l))
			END DO

			! Internal energy !
			CALL EOSEPSILON_NM (primR(irho,j,k,l-1), primR(itau,j,k,l-1), epsR(j,k,l-1))
			CALL EOSEPSILON_NM (primL(irho,j,k,l), primL(itau,j,k,l), epsL(j,k,l))
			CALL EOSSOUNDSPEED (primR(itau,j,k,l-1), primR(irho,j,k,l-1), epsR(j,k,l-1), csR(j,k,l-1))
			CALL EOSSOUNDSPEED (primL(itau,j,k,l), primL(irho,j,k,l), epsL(j,k,l), csL(j,k,l))
	
			! Magnetic field !
			CALL PPMC (wz(l,1:14), bcell(ibx,j,k,l-2), bcell(ibx,j,k,l-1), bcell(ibx,j,k,l), bcell(ibx,j,k,l+1), bcell(ibx,j,k,l+2), primR(ibx,j,k,l-1), primL(ibx,j,k,l))
			CALL PPMC (wz(l,1:14), bcell(iby,j,k,l-2), bcell(iby,j,k,l-1), bcell(iby,j,k,l), bcell(iby,j,k,l+1), bcell(iby,j,k,l+2), primR(iby,j,k,l-1), primL(iby,j,k,l))

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
WRITE(*,*) 'ppmcz = ', REAL(time_end - time_start) / rate
#endif

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate cell average values to interface values using PPM 
! Assumed non-uniform gridding. Using the original PPM algorithm
! See Colella 1984. No steepening or flattening is performed.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PPMC (wgt, vm2, vm1, vc, vp1, vp2, vm_out, vp_out)
!$ACC ROUTINE SEQ 
USE DEFINITION
IMPLICIT NONE

! The input into the subroutine, including conservative variable and input flux function !
REAL*8, INTENT (IN) :: vm2, vm1, vc, vp1, vp2

! The output of the subroutine, the flux at cell boundary !
REAL*8, INTENT (OUT) :: vp_out, vm_out

! Temporal variables !
REAL*8 :: vl, vr
REAL*8 :: deltam1, deltac, deltap1
REAL*8 :: dmp1, dmc, dmm1
REAL*8 :: condition

! Grid size !
REAL*8, DIMENSION(1:14) :: wgt

! Reconstruction coefficients !
REAL*8 :: a0m, a1m, a2m, a3m
REAL*8 :: a0p, a1p, a2p, a3p
REAL*8 :: b0m, b1m
REAL*8 :: b0c, b1c
REAL*8 :: b0p, b1p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign weight !
a0m = wgt(1)
a1m = wgt(2)
a2m = wgt(3)
a3m = wgt(4)
a0p = wgt(5)
a1p = wgt(6)
a2p = wgt(7)
a3p = wgt(8)
b0m = wgt(9)
b1m = wgt(10)
b0c = wgt(11)
b1c = wgt(12)
b0p = wgt(13)
b1p = wgt(14)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Get slopes !
deltap1 = b0p*(vp2 - vp1) + b1p*(vp1 - vc)
deltac = b0c*(vp1 - vc) + b1c*(vc - vm1)
deltam1 = b0m*(vc - vm1) + b1m*(vm1 - vm2)

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
vp_out = vc + a0p*(vp1 - vc) + a1p*(vp1 - vc) + a2p*dmp1 + a3p*dmc
vm_out = vm1 + a0m*(vc - vm1) + a1m*(vc - vm1) + a2m*dmc + a3m*dmm1

! backup !
vl = vm_out
vr = vp_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
