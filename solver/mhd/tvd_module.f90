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

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do reconstruction using the TVD method with the modified van-leer limiter
! along the x-direction, get states at x-interfaces 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE TVDVL_reconx
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
				CALL TVD_VL (x2cen(j-1), x2cen(j), x2cen(j+1), xF2(j-1), xF2(j), dx2(j), & 
										 prim2(i,j-1,k,l), prim2(i,j,k,l), prim2(i,j+1,k,l), & 
										 primR2(i,j-1,k,l), primL2(i,j,k,l))
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
				CALL TVD_VL (x2cen(j-1), x2cen(j), x2cen(j+1), xF2(j-1), xF2(j), dx2(j), & 
									   cs2(j-1,k,l), cs2(j,k,l), cs2(j+1,k,l), & 
										 cs2R(j-1,k,l), cs2L(j,k,l))
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do reconstruction using the TVD method with the modified MC limiter
! along the x-direction, get states at x-interfaces 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE TVDMC_reconx
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
				CALL TVD_MC (x2cen(j-1), x2cen(j), x2cen(j+1), xF2(j-1), xF2(j), dx2(j), & 
										 prim2(i,j-1,k,l), prim2(i,j,k,l), prim2(i,j+1,k,l), & 
										 primR2(i,j-1,k,l), primL2(i,j,k,l))
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
				CALL TVD_MC (x2cen(j-1), x2cen(j), x2cen(j+1), xF2(j-1), xF2(j), dx2(j), & 
									   cs2(j-1,k,l), cs2(j,k,l), cs2(j+1,k,l), & 
										 cs2R(j-1,k,l), cs2L(j,k,l))
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do reconstruction using the TVD method with the modified van-leer (VL) limiter
! along the y-direction, get states at y-interfaces 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE TVDVL_recony
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
				CALL TVD_VL (y2cen(k-1), y2cen(k), y2cen(k+1), yF2(k-1), yF2(k), dy2(k), & 
										 prim2(i,j,k-1,l), prim2(i,j,k,l), prim2(i,j,k+1,l), & 
										 primR2(i,j,k-1,l), primL2(i,j,k,l))
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
				CALL TVD_VL (y2cen(k-1), y2cen(k), y2cen(k+1), yF2(k-1), yF2(k), dy2(k), & 
										 cs2(j,k-1,l), cs2(j,k,l), cs2(j,k+1,l), & 
										 cs2R(j,k-1,l), cs2L(j,k,l))
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do reconstruction using the TVD method with the modified MC limiter
! along the y-direction, get states at y-interfaces 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE TVDMC_recony
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
				CALL TVD_MC (y2cen(k-1), y2cen(k), y2cen(k+1), yF2(k-1), yF2(k), dy2(k), & 
										 prim2(i,j,k-1,l), prim2(i,j,k,l), prim2(i,j,k+1,l), & 
										 primR2(i,j,k-1,l), primL2(i,j,k,l))
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
				CALL TVD_MC (y2cen(k-1), y2cen(k), y2cen(k+1), yF2(k-1), yF2(k), dy2(k), & 
										 cs2(j,k-1,l), cs2(j,k,l), cs2(j,k+1,l), & 
										 cs2R(j,k-1,l), cs2L(j,k,l))
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
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do reconstruction using the TVD method with the modified van leer (VL) limiter
! along the z-direction, get states at z-interfaces 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE TVDVL_reconz
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
				CALL TVD_VL (z2(l-1), z2(l), z2(l+1), zF2(l-1), zF2(l), dz2(l), & 
										 prim2(i,j,k,l-1), prim2(i,j,k,l), prim2(i,j,k,l+1), &
										 primR2(i,j,k,l-1), primL2(i,j,k,l))
			END DO
		END DO
	END DO
END DO
!$OMP END DO

! Dual energy !
IF(dual_energy) THEN
  !$OMP DO COLLAPSE(3) SCHEDULE (STATIC)
	DO l = nz_min_2 - 1, nz_part_2 + 1
		DO k = ny_min_2, ny_part_2
			DO j = nx_min_2, nx_part_2
				eps2R(j,k,l) = primR2(ieps2,j,k,l)/primR2(irho2,j,k,l)
				eps2L(j,k,l) = primL2(ieps2,j,k,l)/primL2(irho2,j,k,l)
				CALL TVD_VL (z2(l-1), z2(l), z2(l+1), zF2(l-1), zF2(l), dz2(l), & 
									   cs2(j,k,l-1), cs2(j,k,l), cs2(j,k,l+1), & 
										 cs2R(j,k,l-1), cs2L(j,k,l))
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do reconstruction using the TVD method with the MC limiter
! along the z-direction, get states at z-interfaces 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE TVDMC_reconz
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
				CALL TVD_MC (z2(l-1), z2(l), z2(l+1), zF2(l-1), zF2(l), dz2(l), & 
										 prim2(i,j,k,l-1), prim2(i,j,k,l), prim2(i,j,k,l+1), &
										 primR2(i,j,k,l-1), primL2(i,j,k,l))
			END DO
		END DO
	END DO
END DO
!$OMP END DO

! Dual energy !
IF(dual_energy) THEN
	!$OMP DO COLLAPSE(3) SCHEDULE (STATIC)
	DO l = nz_min_2 - 1, nz_part_2 + 1
		DO k = ny_min_2, ny_part_2
			DO j = nx_min_2, nx_part_2
				eps2R(j,k,l) = primR2(ieps2,j,k,l)/primR2(irho2,j,k,l)
				eps2L(j,k,l) = primL2(ieps2,j,k,l)/primL2(irho2,j,k,l)
				CALL TVD_MC (z2(l-1), z2(l), z2(l+1), zF2(l-1), zF2(l), dz2(l), & 
									   cs2(j,k,l-1), cs2(j,k,l), cs2(j,k,l+1), & 
										 cs2R(j,k,l-1), cs2L(j,k,l))
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The modified van-leer (VL) slope limiter for getting interface primitive variables 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE TVD_VL (cenm, cenc, cenp, xfm, xfp, dx, vm, vc, vp, vm_out, vp_out)
USE DEFINITION
IMPLICIT NONE

! Input conservative variables !
REAL*8, INTENT (IN) :: cenp, cenm, cenc
REAL*8, INTENT (IN) :: xfp, xfm, dx
REAL*8, INTENT (IN) :: vm, vc, vp

! The reconstructed states at cell boundary !
REAL*8, INTENT (OUT) :: vm_out, vp_out

! Integer !
INTEGER :: i, j, k

! Real !
REAL*8 :: cf, cb, dqm, dq2
REAL*8 :: dqb, dqf
REAL*8 :: slope

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! assign !
dqf = dx*(vp - vc)/(cenp - cenc)
dqb = dx*(vc - vm)/(cenc - cenm)

! Assign !
cf = (cenp - cenc)/(xfp - cenc)
cb = (cenc - cenm)/(cenc - xfm)

! VL limiter !
dq2 = dqf*dqb
dqm = (dq2*(cf*dqb + cb*dqf)/(dqb**2 + dqf**2 + dq2*(cf + cb - 2.0)))
if (dq2 <= 0.0D0) THEN
	dqm = 0.0
END IF

! Output !
vp_out = vc + dqm*(xfp - cenc)/(dx)
vm_out = vc + dqm*(xfm - cenc)/(dx)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The modified MC (monotonized-central) slope limiter for getting interface primitive variables 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE TVD_MC (cenm, cenc, cenp, xfm, xfp, dx, vm, vc, vp, vm_out, vp_out)
USE DEFINITION
IMPLICIT NONE

! Input conservative variables !
REAL*8, INTENT (IN) :: cenp, cenm, cenc
REAL*8, INTENT (IN) :: xfp, xfm, dx
REAL*8, INTENT (IN) :: vm, vc, vp

! The reconstructed states at cell boundary !
REAL*8, INTENT (OUT) :: vm_out, vp_out

! Integer !
INTEGER :: i, j, k

! Real !
REAL*8 :: cf, cb, dqm
REAL*8 :: dqb, dqf, v_tvd
REAL*8 :: slope

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! assign !
dqf = dx*(vp - vc)/(cenp - cenc)
dqb = dx*(vc - vm)/(cenc - cenm)

! Assign !
cf = (cenp - cenc)/(xfp - cenc)
cb = (cenc - cenm)/(cenc - xfm)

! MC limiter !
v_tvd = dqb/dqf
slope = MAX(0.0D0, MIN(0.5D0*(1.0D0 + v_tvd), cf, cb*v_tvd))
dqm = dqf*slope

! Output !
vp_out = vc + dqm*(xfp - cenc)/(dx)
vm_out = vc + dqm*(xfm - cenc)/(dx)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

END MODULE
