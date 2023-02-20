!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This module contains all necessary code for doing the 
! WENO reconstruction. 
! Prototype developed by Wong Ka Wing in 2010 (or before?)
! Merged and systematized by Leung Shing Chi in 2016
! More information about WENO, refer Shu (2000)
!
! This module contains the following subroutines
! 1. subroutine GetConst
! 2. subroutine WENO_r
! 3. subroutine WENO_z
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE WENO_MODULE        
USE DEFINITION
IMPLICIT NONE

! small parameter !
REAL*8, PARAMETER :: smallpara = 1.0D-40

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do reconstruction using WENO interpolation, but along the vertical directions !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE WENO_reconx
USE RIEMANN_MODULE
USE DEFINITION  
IMPLICIT NONE

! Integer !
INTEGER :: j, k, l, i, p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We first interpolate density for NM !
!$OMP PARALLEL
!$OMP DO COLLAPSE(3) SCHEDULE (STATIC)
DO l = nz_min_2, nz_part_2
	DO k = ny_min_2, ny_part_2
		DO j = nx_min_2 - 1, nx_part_2 + 1
			DO i = imin2, imax2
				CALL WENO_MN (dx2(j-2), dx2(j-1), dx2(j), dx2(j+1), dx2(j+2), prim2(i,j-2,k,l), prim2(i,j-1,k,l), prim2(i,j,k,l), prim2(i,j+1,k,l), prim2(i,j+2,k,l), primR2(i,j-1,k,l), primL2(i,j,k,l))
			END DO
		END DO
	END DO
END DO
!$OMP END DO

! Dual energy !
IF(dual_energy) THEN
  !$OMP DO COLLAPSE(3) SCHEDULE (STATIC)
  DO l = nz_min_2, nz_part_2
		DO k = ny_min_2, ny_part_2
			DO j = nx_min_2 - 1, nx_part_2 + 1
				eps2R(j,k,l) = primR2(ieps2,j,k,l)/primR2(irho2,j,k,l)
				eps2L(j,k,l) = primL2(ieps2,j,k,l)/primL2(irho2,j,k,l)
				CALL WENO_MN (dx2(j-2), dx2(j-1), dx2(j), dx2(j+1), dx2(j+2), cs2(j-2,k,l), cs2(j-1,k,l), cs2(j,k,l), cs2(j+1,k,l), cs2(j+2,k,l), cs2R(j-1,k,l), cs2L(j,k,l))
			END DO
		END DO
  END DO
  !$OMP END DO
ELSE
  !$OMP DO COLLAPSE(3) SCHEDULE (STATIC)
  DO l = nz_min_2, nz_part_2
		DO k = ny_min_2, ny_part_2
			DO j = nx_min_2 - 1, nx_part_2
				CALL EOSEPSILON_NM (primR2(irho2,j,k,l), primR2(itau2,j,k,l), eps2R(j,k,l))
				CALL EOSEPSILON_NM (primL2(irho2,j,k,l), primL2(itau2,j,k,l), eps2L(j,k,l))
				CALL EOSSOUNDSPEED (primR2(itau2,j,k,l), primR2(irho2,j,k,l), eps2R(j,k,l), cs2R(j,k,l))
				CALL EOSSOUNDSPEED (primL2(itau2,j,k,l), primL2(irho2,j,k,l), eps2L(j,k,l), cs2L(j,k,l))
			END DO
		END DO
  END DO
  !$OMP END DO
END IF
!$OMP END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do reconstruction using WENO interpolation, but along the vertical directions !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE WENO_recony
USE RIEMANN_MODULE
USE DEFINITION  
IMPLICIT NONE

! Integer !
INTEGER :: j, k, l, i, p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We first interpolate density for NM !
!$OMP PARALLEL
!$OMP DO COLLAPSE(3) SCHEDULE (STATIC)
DO l = nz_min_2, nz_part_2
	DO k = ny_min_2 - 1, ny_part_2 + 1
		DO j = nx_min_2, nx_part_2
			DO i = imin2, imax2
				CALL WENO_MN (dy2(k-2), dy2(k-1), dy2(k), dy2(k+1), dy2(k+2), prim2(i,j,k-2,l), prim2(i,j,k-1,l), prim2(i,j,k,l), prim2(i,j,k+1,l), prim2(i,j,k+2,l), primR2(i,j,k-1,l), primL2(i,j,k,l))
			END DO
		END DO
	END DO
END DO
!$OMP END DO

! Dual energy !
IF(dual_energy) THEN
  !$OMP DO COLLAPSE(3) SCHEDULE (STATIC)
  DO l = nz_min_2, nz_part_2
		DO k = ny_min_2 - 1, ny_part_2 + 1
			DO j = nx_min_2, nx_part_2
				eps2R(j,k,l) = primR2(ieps2,j,k,l)/primR2(irho2,j,k,l)
				eps2L(j,k,l) = primL2(ieps2,j,k,l)/primL2(irho2,j,k,l)
				CALL WENO_MN (dy2(k-2), dy2(k-1), dy2(k), dy2(k+1), dy2(k+2), cs2(j,k-2,l), cs2(j,k-1,l), cs2(j,k,l), cs2(j,k+1,l), cs2(j,k+2,l), cs2R(j,k-1,l), cs2L(j,k,l))
			END DO
		END DO
  END DO
  !$OMP END DO
ELSE
  !$OMP DO COLLAPSE(3) SCHEDULE (STATIC)
  DO l = nz_min_2, nz_part_2
		DO k = ny_min_2 - 1, ny_part_2
			DO j = nx_min_2, nx_part_2
				CALL EOSEPSILON_NM (primR2(irho2,j,k,l), primR2(itau2,j,k,l), eps2R(j,k,l))
				CALL EOSEPSILON_NM (primL2(irho2,j,k,l), primL2(itau2,j,k,l), eps2L(j,k,l))
				CALL EOSSOUNDSPEED (primR2(itau2,j,k,l), primR2(irho2,j,k,l), eps2R(j,k,l), cs2R(j,k,l))
				CALL EOSSOUNDSPEED (primL2(itau2,j,k,l), primL2(irho2,j,k,l), eps2L(j,k,l), cs2L(j,k,l))
			END DO
		END DO
  END DO
  !$OMP END DO
END IF
!$OMP END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do reconstruction using WENO interpolation, but along the vertical directions !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE WENO_reconz
USE RIEMANN_MODULE
USE DEFINITION  
IMPLICIT NONE

! Integer !
INTEGER :: j, k, l, i, p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We first interpolate density for NM !
!$OMP PARALLEL
!$OMP DO COLLAPSE(3) SCHEDULE (STATIC)
DO l = nz_min_2 - 1, nz_part_2 + 1
	DO k = ny_min_2, ny_part_2
		DO j = nx_min_2, nx_part_2
			DO i = imin2, imax2
				CALL WENO_MN (dz2(l-2), dz2(l-1), dz2(l), dz2(l+1), dz2(l+2), prim2(i,j,k,l-2), prim2(i,j,k,l-1), prim2(i,j,k,l), prim2(i,j,k,l+1), prim2(i,j,k,l+2), primR2(i,j,k,l-1), primL2(i,j,k,l))
			END DO
		END DO
	END DO
END DO
!$OMP END DO

! Dual energy !
IF(dual_energy) THEN
  !$OMP DO COLLAPSE(3) SCHEDULE (STATIC)
  DO l = nz_min_2 - 1, nz_part_2 + 1
		DO k = ny_min_2 , ny_part_2
			DO j = nx_min_2, nx_part_2
				eps2R(j,k,l) = primR2(ieps2,j,k,l)/primR2(irho2,j,k,l)
				eps2L(j,k,l) = primL2(ieps2,j,k,l)/primL2(irho2,j,k,l)
				CALL WENO_MN (dz2(l-2), dz2(l-1), dz2(l), dz2(l+1), dz2(l+2), cs2(j,k,l-2), cs2(j,k,l-1), cs2(j,k,l), cs2(j,k,l+1), cs2(j,k,l+2), cs2R(j,k,l-1), cs2L(j,k,l))
			END DO
		END DO
  END DO
  !$OMP END DO
ELSE
  !$OMP DO COLLAPSE(3) SCHEDULE (STATIC)
  DO l = nz_min_2 - 1, nz_part_2
		DO k = ny_min_2, ny_part_2
			DO j = nx_min_2, nx_part_2
				CALL EOSEPSILON_NM (primR2(irho2,j,k,l), primR2(itau2,j,k,l), eps2R(j,k,l))
				CALL EOSEPSILON_NM (primL2(irho2,j,k,l), primL2(itau2,j,k,l), eps2L(j,k,l))
				CALL EOSSOUNDSPEED (primR2(itau2,j,k,l), primR2(irho2,j,k,l), eps2R(j,k,l), cs2R(j,k,l))
				CALL EOSSOUNDSPEED (primL2(itau2,j,k,l), primL2(irho2,j,k,l), eps2L(j,k,l), cs2L(j,k,l))
			END DO
		END DO
  END DO
  !$OMP END DO
END IF
!$OMP END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WENO-Z reconstruction of primitive variables on cell interfaces, assuming 		    !
! non-equal mesh size For details, please refer to the textbook with ISBN 3-540-65893-9 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE WENO_MN (dm2, dm1, dc, dp1, dp2, vm2, vm1, vc, vp1, vp2, vm_out, vp_out)
USE DEFINITION
IMPLICIT NONE

! Input conservative variables !
REAL*8, INTENT (IN) :: dm2, dm1, dc, dp1, dp2
REAL*8, INTENT (IN) :: vm2, vm1, vc, vp1, vp2

! The reconstructed states at cell boundary !
REAL*8, INTENT (OUT) :: vm_out, vp_out

! Integer !
INTEGER :: i, j, k

! Tempeorary parameter !
REAL*8 :: tau5, temp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Patches from FEST-3D !
REAL*8 :: P_weno_1, P_weno_2, P_weno_3 
REAL*8 :: B_weno_1, B_weno_2, B_weno_3
REAL*8 :: w_weno_1, w_weno_2, w_weno_3
REAL*8 :: g_weno_1, g_weno_2, g_weno_3

! Real variables !
REAL*8 :: U11_weno
REAL*8 :: U00_weno
REAL*8 :: U21_weno
REAL*8 :: U10_weno
REAL*8 :: U01_weno
REAL*8 :: U12_weno
REAL*8 :: alpha12_weno
REAL*8 :: alpha01_weno
REAL*8 :: alpha10_weno
REAL*8 :: alpha21_weno

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Linear weight !
g_weno_1 = 1.0D0/10.0D0
g_weno_2 = 6.0D0/10.0D0
g_weno_3 = 3.0D0/10.0D0

! Assign alpha !
alpha12_weno = dp2/(dp1 + dp2)
alpha01_weno = dp1/(dc + dp1)
alpha10_weno = dc/(dm1 + dc)
alpha21_weno = dm1/(dm2 + dm1)

! Assign slopes !
U01_weno = (1.0D0 - alpha01_weno)*vc + alpha01_weno*vp1
U12_weno = (1.0D0 - alpha12_weno)*vp1 + alpha12_weno*vp2
U10_weno = (1.0D0 - alpha10_weno)*vm1 + alpha10_weno*vc
U21_weno = (1.0D0 - alpha21_weno)*vm2 + alpha21_weno*vm1
U00_weno = vm1 + (1.0D0 - alpha21_weno)*(vm1 - vm2)
U11_weno = vp1 + alpha12_weno*(vp1 - vp2)

! Assign polynominals !
P_weno_1 = ( 6.0D0*vc - 1.0D0*U10_weno - 2.0D0*U00_weno)/3.0D0
P_weno_2 = (- 1.0D0*U10_weno + 2.0D0*vc + 2.0D0*U01_weno)/3.0D0
P_weno_3 = ( 2.0D0*U01_weno + 2.0D0*vp1 - 1.0D0*U12_weno)/3.0D0

! Assign smoothness indicator !
B_weno_1 = (13.0D0/12.0D0)*(2.0D0*U10_weno - 2.0*U00_weno)**2 + & 
	    (1.0D0/4.0D0)*(4.0D0*vc - 2.0D0*U10_weno - 2.0D0*U00_weno)**2
B_weno_2 = (13.0D0/12.0D0)*(2.0D0*U10_weno - 4.0D0*vc + 2.0D0*U01_weno)**2 + &
	    (1.0D0/4.0D0)*(2.0D0*U01_weno - 2.0D0*U10_weno)**2
B_weno_3 = (13.0D0/12.0D0)*(2.0D0*U01_weno - 4.0D0*vp1 + 2.0D0*U12_weno)**2 + & 
	    (1.0D0/4.0D0)*(8.0D0*vp1 - 6.0D0*U01_weno - 2.0D0*U12_weno)**2

! WENO-Z !
tau5 = abs(B_weno_1 - B_weno_3)

! Weight !
w_weno_1 = g_weno_1*(1.0D0 + (tau5/(smallpara + B_weno_1))**2)
w_weno_2 = g_weno_2*(1.0D0 + (tau5/(smallpara + B_weno_2))**2)
w_weno_3 = g_weno_3*(1.0D0 + (tau5/(smallpara + B_weno_3))**2)

! Output !
vp_out = (w_weno_1*P_weno_1 + w_weno_2*P_weno_2 + w_weno_3*P_weno_3)/(w_weno_1 + w_weno_2 + w_weno_3)

! Do the same for the right states !
P_weno_1 = (6.0D0*vc - 1.0D0*U01_weno - 2.0D0*U11_weno)/3.0D0
P_weno_2 = (-1.0D0*U01_weno + 2.0D0*vc + 2.0D0*U10_weno)/3.0D0
P_weno_3 = (2.0D0*U10_weno + 2.0D0*vm1 - 1.0D0*U21_weno)/3.0D0

! Assign smoothness indicator !
B_weno_1 = (13.0D0/12.0D0)*(2.0D0*U01_weno - 2.0D0*U11_weno)**2 + & 
	    (1.0D0/4.0D0)*(4.0D0*vc - 2.0D0*U01_weno - 2.0D0*U11_weno)**2
B_weno_2 = (13.0D0/12.0D0)*(2.0D0*U01_weno - 4.0D0*vc + 2.0D0*U10_weno)**2 + &
	    (1.0D0/4.0D0)*(-2.0D0*u01_weno + 2.0D0*U10_weno)**2
B_weno_3 = (13.0D0/12.0D0)*(2.0D0*U10_weno - 4.0D0*vm1 + 2.0D0*U21_weno)**2 + & 
	    (1.0D0/4.0D0)*(-6.0D0*U10_weno + 8.0D0*vm1 - 2.0D0*U21_weno)**2

! WENO-Z !
tau5 = abs(B_weno_1 - B_weno_3)

! Weight !
w_weno_1 = g_weno_1*(1.0D0 + (tau5/(smallpara + B_weno_1))**2)
w_weno_2 = g_weno_2*(1.0D0 + (tau5/(smallpara + B_weno_2))**2)
w_weno_3 = g_weno_3*(1.0D0 + (tau5/(smallpara + B_weno_3))**2)

! Output !
vm_out = (w_weno_1*P_weno_1 + w_weno_2*P_weno_2 + w_weno_3*P_weno_3)/(w_weno_1 + w_weno_2 + w_weno_3)

END SUBROUTINE

END MODULE
