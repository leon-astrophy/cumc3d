!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This module contains all necessary code for doing the 
! WENO reconstruction. 
! Prototype developed by Wong Ka Wing in 2010 (or before?)
! Merged and systematized by Leung Shing Chi in 2016
! More information about WENO, refer Shu (2000)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE WENO_MODULE        
USE DEFINITION
IMPLICIT NONE

! small parameter !
REAL*8, PARAMETER :: smallpara = 1.0D-40

! Linear weight !
REAL*8, PARAMETER :: g_weno_1 = 1.0D0/10.0D0
REAL*8, PARAMETER :: g_weno_2 = 6.0D0/10.0D0
REAL*8, PARAMETER :: g_weno_3 = 3.0D0/10.0D0

! Reconstruction coefficents !
REAL*8, PARAMETER :: r_weno_1 = 13.0D0/12.0D0
REAL*8, PARAMETER :: r_weno_2 = 1.0D0/4.0D0

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Setup WENO interpolation weight 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE WENO_WEIGHT
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k, l

! Real !
REAL*8 :: dm2, dm1, dc, dp1, dp2
REAL*8 :: alpha12_weno, alpha01_weno
REAL*8 :: alpha10_weno, alpha21_weno

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set up for x, y, and z directional reconstruction !

ALLOCATE(wx(-2:nx+3,1:10))
ALLOCATE(wy(-2:ny+3,1:10))
ALLOCATE(wz(-2:nz+3,1:10))

! Assign alpha !
alpha12_weno = dp2/(dp1 + dp2)
alpha01_weno = dp1/(dc + dp1)
alpha10_weno = dc/(dm1 + dc)
alpha21_weno = dm1/(dm2 + dm1)

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
	CALL WENO_INTERPOLANT(dm2, dm1, dc, dp1, dp2, wx(j,1:10))

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
	CALL WENO_INTERPOLANT(dm2, dm1, dc, dp1, dp2, wy(j,1:10))

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
	CALL WENO_INTERPOLANT(dm2, dm1, dc, dp1, dp2, wz(j,1:10))

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Get WENO interpolant weight from the uneven grid size dx
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE WENO_INTERPOLANT(dm2, dm1, dc, dp1, dp2, wgt_out)
USE DEFINITION
IMPLICIT NONE

! Input !
REAL*8, INTENT(IN) :: dm2, dm1, dc, dp1, dp2
REAL*8, DIMENSION(1:10) :: wgt_out

! Real !
REAL*8 :: alpha12_weno 
REAL*8 :: alpha01_weno
REAL*8 :: alpha10_weno 
REAL*8 :: alpha21_weno 

! Assign alpha !
alpha12_weno = dp2/(dp1 + dp2)
alpha01_weno = dp1/(dc + dp1)
alpha10_weno = dc/(dm1 + dc)
alpha21_weno = dm1/(dm2 + dm1)

! Assign weight !
wgt_out(1) = (1.0D0 - alpha01_weno)
wgt_out(2) = alpha01_weno
wgt_out(3) = (1.0D0 - alpha12_weno)
wgt_out(4) = alpha12_weno
wgt_out(5) = (1.0D0 - alpha10_weno)
wgt_out(6) = alpha10_weno
wgt_out(7) = (1.0D0 - alpha21_weno)
wgt_out(8) = alpha21_weno
wgt_out(9) = (1.0D0 - alpha21_weno)
wgt_out(10) = alpha12_weno

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do reconstruction using WENO interpolation, but along the vertical directions !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE WENO_reconx
!$ACC ROUTINE (EOSEPSILON_NM) SEQ
!$ACC ROUTINE (EOSSOUNDSPEED) SEQ
!$ACC ROUTINE (WENO) SEQ
USE RIEMANN_MODULE
USE MHD_MODULE
USE DEFINITION  
IMPLICIT NONE

! Integer !
INTEGER :: j, k, l, i, p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We first interpolate primitive variables for NM !
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE (STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT (PRESENT)
DO l = 1, nz
	DO k = 1, ny
		DO j = 0, nx + 1 

			! NM hydro variables !
			DO i = imin, ibx - 1
				CALL WENO (wx(j,1:10), prim(i,j-2,k,l), prim(i,j-1,k,l), prim(i,j,k,l), prim(i,j+1,k,l), prim(i,j+2,k,l), primR(i,j-1,k,l), primL(i,j,k,l))
			END DO

			! Internal energy !
			CALL EOSEPSILON_NM (primR(irho,j-1,k,l), primR(itau,j-1,k,l), epsR(j-1,k,l))
			CALL EOSEPSILON_NM (primL(irho,j,k,l), primL(itau,j,k,l), epsL(j,k,l))
			CALL EOSSOUNDSPEED (primR(itau,j-1,k,l), primR(irho,j-1,k,l), epsR(j-1,k,l), csR(j-1,k,l))
			CALL EOSSOUNDSPEED (primL(itau,j,k,l), primL(irho,j,k,l), epsL(j,k,l), csL(j,k,l))

			! magnetic field !
			CALL WENO (wx(j,1:10), bcell(iby,j-2,k,l), bcell(iby,j-1,k,l), bcell(iby,j,k,l), bcell(iby,j+1,k,l), bcell(iby,j+2,k,l), primR(iby,j-1,k,l), primL(iby,j,k,l))
			CALL WENO (wx(j,1:10), bcell(ibz,j-2,k,l), bcell(ibz,j-1,k,l), bcell(ibz,j,k,l), bcell(ibz,j+1,k,l), bcell(ibz,j+2,k,l), primR(ibz,j-1,k,l), primL(ibz,j,k,l))

			! Special treatment for normal field 
			primR(ibx,j-1,k,l) = prim(ibx,j-1,k,l)
			primL(ibx,j,k,l) = prim(ibx,j,k,l)

		END DO
	END DO
END DO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do reconstruction using WENO interpolation, but along the vertical directions !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE WENO_recony
!$ACC ROUTINE (EOSEPSILON_NM) SEQ
!$ACC ROUTINE (EOSSOUNDSPEED) SEQ
!$ACC ROUTINE (WENO) SEQ
USE RIEMANN_MODULE
USE MHD_MODULE
USE DEFINITION  
IMPLICIT NONE

! Integer !
INTEGER :: j, k, l, i, p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We first interpolate primitive variables for NM !
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE (STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT (PRESENT)
DO l = 1, nz
	DO k = 0, ny + 1 
		DO j = 1, nx

			! NM hydro variables !
			DO i = imin, ibx - 1
				CALL WENO (wy(k,1:10), prim(i,j,k-2,l), prim(i,j,k-1,l), prim(i,j,k,l), prim(i,j,k+1,l), prim(i,j,k+2,l), primR(i,j,k-1,l), primL(i,j,k,l))
			END DO

			! Internal energy !
			CALL EOSEPSILON_NM (primR(irho,j,k-1,l), primR(itau,j,k-1,l), epsR(j,k-1,l))
			CALL EOSEPSILON_NM (primL(irho,j,k,l), primL(itau,j,k,l), epsL(j,k,l))
			CALL EOSSOUNDSPEED (primR(itau,j,k-1,l), primR(irho,j,k-1,l), epsR(j,k-1,l), csR(j,k-1,l))
			CALL EOSSOUNDSPEED (primL(itau,j,k,l), primL(irho,j,k,l), epsL(j,k,l), csL(j,k,l))

			! Magnetic fields !
			CALL WENO (wy(k,1:10), bcell(ibx,j,k-2,l), bcell(ibx,j,k-1,l), bcell(ibx,j,k,l), bcell(ibx,j,k+1,l), bcell(ibx,j,k+2,l), primR(ibx,j,k-1,l), primL(ibx,j,k,l))
			CALL WENO (wy(k,1:10), bcell(ibz,j,k-2,l), bcell(ibz,j,k-1,l), bcell(ibz,j,k,l), bcell(ibz,j,k+1,l), bcell(ibz,j,k+2,l), primR(ibz,j,k-1,l), primL(ibz,j,k,l))

			! special treatment for formal fields !
			primR(iby,j,k-1,l) = prim(iby,j,k-1,l)
			primL(iby,j,k,l) = prim(iby,j,k,l)	
										 
		END DO
	END DO
END DO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do reconstruction using WENO interpolation, but along the vertical directions !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE WENO_reconz
!$ACC ROUTINE (EOSEPSILON_NM) SEQ
!$ACC ROUTINE (EOSSOUNDSPEED) SEQ
!$ACC ROUTINE (WENO) SEQ
USE RIEMANN_MODULE
USE MHD_MODULE
USE DEFINITION  
IMPLICIT NONE

! Integer !
INTEGER :: j, k, l, i, p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We first interpolate primitive variables for NM !
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE (STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT (PRESENT)
DO l = 0, nz + 1 
	DO k = 1, ny
		DO j = 1, nx

			! NM hydro variables !
			DO i = imin, ibx - 1
				CALL WENO (wz(l,1:10), prim(i,j,k,l-2), prim(i,j,k,l-1), prim(i,j,k,l), prim(i,j,k,l+1), prim(i,j,k,l+2), primR(i,j,k,l-1), primL(i,j,k,l))
			END DO

			! Internal energy !
			CALL EOSEPSILON_NM (primR(irho,j,k,l-1), primR(itau,j,k,l-1), epsR(j,k,l-1))
			CALL EOSEPSILON_NM (primL(irho,j,k,l), primL(itau,j,k,l), epsL(j,k,l))
			CALL EOSSOUNDSPEED (primR(itau,j,k,l-1), primR(irho,j,k,l-1), epsR(j,k,l-1), csR(j,k,l-1))
			CALL EOSSOUNDSPEED (primL(itau,j,k,l), primL(irho,j,k,l), epsL(j,k,l), csL(j,k,l))
	
			! Magnetic field !
			CALL WENO (wz(l,1:10), bcell(ibx,j,k,l-2), bcell(ibx,j,k,l-1), bcell(ibx,j,k,l), bcell(ibx,j,k,l+1), bcell(ibx,j,k,l+2), primR(ibx,j,k,l-1), primL(ibx,j,k,l))
			CALL WENO (wz(l,1:10), bcell(iby,j,k,l-2), bcell(iby,j,k,l-1), bcell(iby,j,k,l), bcell(iby,j,k,l+1), bcell(iby,j,k,l+2), primR(iby,j,k,l-1), primL(iby,j,k,l))

			! special treatment for normal field !
			primR(ibz,j,k,l-1) = prim(ibz,j,k,l-1)
			primL(ibz,j,k,l) = prim(ibz,j,k,l)										

		END DO
	END DO
END DO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! WENO-Z reconstruction of primitive variables on cell interfaces, assuming uneven mesh
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE WENO (wgt, vm2, vm1, vc, vp1, vp2, vm_out, vp_out)
!$ACC ROUTINE SEQ 
USE DEFINITION
IMPLICIT NONE

! Input conservative variables !
REAL*8, INTENT (IN) :: vm2, vm1, vc, vp1, vp2

! The reconstructed states at cell boundary !
REAL*8, INTENT (OUT) :: vm_out, vp_out

! Weight !
REAL*8, DIMENSION(1:10) :: wgt

! Integer !
INTEGER :: i, j, k

! Tempeorary parameter !
REAL*8 :: tau5, temp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Patches from FEST-3D !
REAL*8 :: P_weno_1, P_weno_2, P_weno_3 
REAL*8 :: B_weno_1, B_weno_2, B_weno_3
REAL*8 :: w_weno_1, w_weno_2, w_weno_3

! Real variables !
REAL*8 :: U11_weno
REAL*8 :: U00_weno
REAL*8 :: U21_weno
REAL*8 :: U10_weno
REAL*8 :: U01_weno
REAL*8 :: U12_weno

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign slopes !
U01_weno = wgt(1)*vc + wgt(2)*vp1
U12_weno = wgt(3)*vp1 + wgt(4)*vp2
U10_weno = wgt(5)*vm1 + wgt(6)*vc
U21_weno = wgt(7)*vm2 + wgt(8)*vm1
U00_weno = vm1 + wgt(9)*(vm1 - vm2)
U11_weno = vp1 + wgt(10)*(vp1 - vp2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do the same for the right states !
P_weno_1 = (6.0D0*vc - 1.0D0*U10_weno - 2.0D0*U00_weno)/3.0D0
P_weno_2 = (-1.0D0*U10_weno + 2.0D0*vc + 2.0D0*U01_weno)/3.0D0
P_weno_3 = (2.0D0*U01_weno + 2.0D0*vp1 - 1.0D0*U12_weno)/3.0D0

! Assign smoothness indicator !
B_weno_1 = r_weno_1*(2.0D0*U10_weno - 2.0D0*U00_weno)**2 + r_weno_2*(4.0D0*vc - 2.0D0*U10_weno - 2.0D0*U00_weno)**2
B_weno_2 = r_weno_1*(2.0D0*U10_weno - 4.0D0*vc + 2.0D0*U01_weno)**2 + r_weno_2*(-2.0D0*U10_weno + 2.0D0*U01_weno)**2
B_weno_3 = r_weno_1*(2.0D0*U01_weno - 4.0D0*vp1 + 2.0D0*U12_weno)**2 + r_weno_2*(-6.0D0*U01_weno + 8.0D0*vp1 - 2.0D0*U12_weno)**2

! WENO-Z !
tau5 = abs(B_weno_1 - B_weno_3)

! Weight !
w_weno_1 = g_weno_1*(1.0D0 + (tau5/(smallpara + B_weno_1))**2)
w_weno_2 = g_weno_2*(1.0D0 + (tau5/(smallpara + B_weno_2))**2)
w_weno_3 = g_weno_3*(1.0D0 + (tau5/(smallpara + B_weno_3))**2)

! Output !
vp_out = (w_weno_1*P_weno_1 + w_weno_2*P_weno_2 + w_weno_3*P_weno_3)/(w_weno_1 + w_weno_2 + w_weno_3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign polynominals !
P_weno_1 = ( 6.0D0*vc - 1.0D0*U01_weno - 2.0D0*U11_weno)/3.0D0
P_weno_2 = (- 1.0D0*U01_weno + 2.0D0*vc + 2.0D0*U10_weno)/3.0D0
P_weno_3 = ( 2.0D0*U10_weno + 2.0D0*vm1 - 1.0D0*U21_weno)/3.0D0

! Assign smoothness indicator !
B_weno_1 = r_weno_1*(2.0D0*U01_weno - 2.0*U11_weno)**2 + r_weno_2*(4.0D0*vc - 2.0D0*U01_weno - 2.0D0*U11_weno)**2
B_weno_2 = r_weno_1*(2.0D0*U01_weno - 4.0D0*vc + 2.0D0*U10_weno)**2 + r_weno_2*(2.0D0*U01_weno - 2.0D0*U10_weno)**2
B_weno_3 = r_weno_1*(2.0D0*U10_weno - 4.0D0*vm1 + 2.0D0*U21_weno)**2 + r_weno_2*(8.0D0*vm1 - 6.0D0*U10_weno - 2.0D0*U21_weno)**2

! WENO-Z !
tau5 = abs(B_weno_1 - B_weno_3)

! Weight !
w_weno_1 = g_weno_1*(1.0D0 + (tau5/(smallpara + B_weno_1))**2)
w_weno_2 = g_weno_2*(1.0D0 + (tau5/(smallpara + B_weno_2))**2)
w_weno_3 = g_weno_3*(1.0D0 + (tau5/(smallpara + B_weno_3))**2)

! Output !
vm_out = (w_weno_1*P_weno_1 + w_weno_2*P_weno_2 + w_weno_3*P_weno_3)/(w_weno_1 + w_weno_2 + w_weno_3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE
