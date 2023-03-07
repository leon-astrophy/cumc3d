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
				CALL TVD_VL (dx2(j-1), dx2(j), dx2(j+1), prim2(i,j-1,k,l), prim2(i,j,k,l), prim2(i,j+1,k,l), primR2(i,j-1,k,l), primL2(i,j,k,l))
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
				CALL TVD_VL (dx2(j-1), dx2(j), dx2(j+1), cs2(j-1,k,l), cs2(j,k,l), cs2(j+1,k,l), cs2R(j-1,k,l), cs2L(j,k,l))
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
				CALL TVD_MC (dx2(j-1), dx2(j), dx2(j+1), prim2(i,j-1,k,l), prim2(i,j,k,l), prim2(i,j+1,k,l), primR2(i,j-1,k,l), primL2(i,j,k,l))
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
				CALL TVD_MC (dx2(j-1), dx2(j), dx2(j+1), cs2(j-1,k,l), cs2(j,k,l), cs2(j+1,k,l), cs2R(j-1,k,l), cs2L(j,k,l))
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
				CALL TVD_VL (dy2(k-1), dy2(k), dy2(k+1), prim2(i,j,k-1,l), prim2(i,j,k,l), prim2(i,j,k+1,l), primR2(i,j,k-1,l), primL2(i,j,k,l))
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
				CALL TVD_VL (dy2(k-1), dy2(k), dy2(k+1), cs2(j,k-1,l), cs2(j,k,l), cs2(j,k+1,l), cs2R(j,k-1,l), cs2L(j,k,l))
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
				CALL TVD_MC (dy2(k-1), dy2(k), dy2(k+1), prim2(i,j,k-1,l), prim2(i,j,k,l), prim2(i,j,k+1,l), primR2(i,j,k-1,l), primL2(i,j,k,l))
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
				CALL TVD_MC (dy2(k-1), dy2(k), dy2(k+1), cs2(j,k-1,l), cs2(j,k,l), cs2(j,k+1,l), cs2R(j,k-1,l), cs2L(j,k,l))
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
				CALL TVD_VL (dz2(l-1), dz2(l), dz2(l+1), prim2(i,j,k,l-1), prim2(i,j,k,l), prim2(i,j,k,l+1), primR2(i,j,k,l-1), primL2(i,j,k,l))
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
				CALL TVD_VL (dz2(l-1), dz2(l), dz2(l+1), cs2(j,k,l-1), cs2(j,k,l), cs2(j,k,l+1), cs2R(j,k,l-1), cs2L(j,k,l))
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
				CALL TVD_MC (dz2(l-1), dz2(l), dz2(l+1), prim2(i,j,k,l-1), prim2(i,j,k,l), prim2(i,j,k,l+1), primR2(i,j,k,l-1), primL2(i,j,k,l))
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
				CALL TVD_MC (dz2(l-1), dz2(l), dz2(l+1), cs2(j,k,l-1), cs2(j,k,l), cs2(j,k,l+1), cs2R(j,k,l-1), cs2L(j,k,l))
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
SUBROUTINE TVD_VL (dm1, dc, dp1, vm1, vc, vp1, vm_out, vp_out)
USE DEFINITION
IMPLICIT NONE

! Input conservative variables !
REAL*8, INTENT (IN) :: dm1, dc, dp1
REAL*8, INTENT (IN) :: vm1, vc, vp1

! The reconstructed states at cell boundary !
REAL*8, INTENT (OUT) :: vm_out, vp_out

! Integer !
INTEGER :: i, j, k

! Real !
REAL*8 :: dql, dqr, dq2, qc
REAL*8 :: cf, cb, dqm
REAL*8 :: dqb, dqf, v_tvd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign !
dql = vc - vm1
dqr = vp1 - vc
qc = vc

! assign !
dqf = dqr*dc/(0.5D0*(dp1 + dc))
dqb = dql*dc/(0.5D0*(dc + dm1))
dq2 = dqf*dqb

! Assign !
cf = (dp1 + dc)/(dc)
cb = (dc + dm1)/(dc)

! Slope !
dqm = (dq2*(cf*dqb + cb*dqf)/(dqb**2 + dqf**2 + dq2*(cf + cb - 2.0)))
if (dq2 <= 0.0D0) THEN
	dqm = 0.0
END IF

! MC limiter !
!v_tvd = dqb/dqf
!dqm = dqf*MAX(0.0D0, MIN(0.5D0*(1.0D0 + v_tvd),  cf, cb*v_tvd))

! Output !
vp_out = vc + 0.5D0*dqm
vm_out = vc - 0.5D0*dqm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The modified MC (monotonized-central) slope limiter for getting interface primitive variables 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE TVD_MC (dm1, dc, dp1, vm1, vc, vp1, vm_out, vp_out)
USE DEFINITION
IMPLICIT NONE

! Input conservative variables !
REAL*8, INTENT (IN) :: dm1, dc, dp1
REAL*8, INTENT (IN) :: vm1, vc, vp1

! The reconstructed states at cell boundary !
REAL*8, INTENT (OUT) :: vm_out, vp_out

! Integer !
INTEGER :: i, j, k

! Real !
REAL*8 :: dql, dqr, dq2, qc
REAL*8 :: cf, cb, dqm
REAL*8 :: dqb, dqf, v_tvd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign !
dql = vc - vm1
dqr = vp1 - vc
qc = vc

! assign !
dqf = dqr*dc/(0.5D0*(dp1 + dc))
dqb = dql*dc/(0.5D0*(dc + dm1))
dq2 = dqf*dqb

! Assign !
cf = (dp1 + dc)/(dc)
cb = (dc + dm1)/(dc)

! MC limiter !
v_tvd = dqb/dqf
dqm = dqf*MAX(0.0D0, MIN(0.5D0*(1.0D0 + v_tvd),  cf, cb*v_tvd))

! Output !
vp_out = vc + 0.5D0*dqm
vm_out = vc - 0.5D0*dqm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

END MODULE
