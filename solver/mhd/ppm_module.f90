!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the PPM (Piecewise Parabolic Method) Module that reconstruct   
! interface values of primitive variables. Three different appoarch 
! can be chosen, and they differ by how to handle monotone conditions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE PPM_MODULE
USE DEFINITION
IMPLICIT NONE

contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do reconstruction using PPM interpolation, with the original Colella algorithm
! Reconstruct along the x-direction, and get states at the x-interfaces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PPM1984_reconx
USE MHD_MODULE
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
				CALL PPM_1984 (dx2(j-2), dx2(j-1), dx2(j), dx2(j+1), dx2(j+2), prim2(i,j-2,k,l), prim2(i,j-1,k,l), prim2(i,j,k,l), prim2(i,j+1,k,l), prim2(i,j+2,k,l), primR2(i,j-1,k,l), primL2(i,j,k,l))
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
				CALL PPM_1984 (dx2(j-2), dx2(j-1), dx2(j), dx2(j+1), dx2(j+2), cs2(j-2,k,l), cs2(j-1,k,l), cs2(j,k,l), cs2(j+1,k,l), cs2(j+2,k,l), cs2R(j-1,k,l), cs2L(j,k,l))
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
! Do reconstruction using PPM interpolation, with the Mignone 2011 algorithm
! Reconstruct along the x-direction, and get states at the x-interfaces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PPM2011_reconx
USE MHD_MODULE
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
				CALL PPM_2011 (dx2(j-2), dx2(j-1), dx2(j), dx2(j+1), dx2(j+2), prim2(i,j-2,k,l), prim2(i,j-1,k,l), prim2(i,j,k,l), prim2(i,j+1,k,l), prim2(i,j+2,k,l), primR2(i,j-1,k,l), primL2(i,j,k,l))
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
				CALL PPM_2011 (dx2(j-2), dx2(j-1), dx2(j), dx2(j+1), dx2(j+2), cs2(j-2,k,l), cs2(j-1,k,l), cs2(j,k,l), cs2(j+1,k,l), cs2(j+2,k,l), cs2R(j-1,k,l), cs2L(j,k,l))
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
! Do reconstruction using PPM interpolation, with the Mignone 2014 algorithm
! Reconstruct along the x-direction, and get states at the x-interfaces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PPM2014_reconx
USE MHD_MODULE
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
				CALL PPM_2014 (dx2(j-2), dx2(j-1), dx2(j), dx2(j+1), dx2(j+2), prim2(i,j-2,k,l), prim2(i,j-1,k,l), prim2(i,j,k,l), prim2(i,j+1,k,l), prim2(i,j+2,k,l), primR2(i,j-1,k,l), primL2(i,j,k,l))
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
				CALL PPM_2014 (dx2(j-2), dx2(j-1), dx2(j), dx2(j+1), dx2(j+2), cs2(j-2,k,l), cs2(j-1,k,l), cs2(j,k,l), cs2(j+1,k,l), cs2(j+2,k,l), cs2R(j-1,k,l), cs2L(j,k,l))
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
! Do reconstruction using PPM interpolation, with the original Colella algorithm
! Reconstruct along the y-direction, and get states at the y-interfaces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PPM1984_recony
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
				CALL PPM_1984 (dy2(k-2), dy2(k-1), dy2(k), dy2(k+1), dy2(k+2), prim2(i,j,k-2,l), prim2(i,j,k-1,l), prim2(i,j,k,l), prim2(i,j,k+1,l), prim2(i,j,k+2,l), primR2(i,j,k-1,l), primL2(i,j,k,l))
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
				CALL PPM_1984 (dy2(k-2), dy2(k-1), dy2(k), dy2(k+1), dy2(k+2), cs2(j,k-2,l), cs2(j,k-1,l), cs2(j,k,l), cs2(j,k+1,l), cs2(j,k+2,l), cs2R(j,k-1,l), cs2L(j,k,l))
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
! Do reconstruction using PPM interpolation, with the Mignone 2011 algorithm
! Reconstruct along the y-direction, and get states at the y-interfaces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PPM2011_recony
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
				CALL PPM_2011 (dy2(k-2), dy2(k-1), dy2(k), dy2(k+1), dy2(k+2), prim2(i,j,k-2,l), prim2(i,j,k-1,l), prim2(i,j,k,l), prim2(i,j,k+1,l), prim2(i,j,k+2,l), primR2(i,j,k-1,l), primL2(i,j,k,l))
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
				CALL PPM_2011 (dy2(k-2), dy2(k-1), dy2(k), dy2(k+1), dy2(k+2), cs2(j,k-2,l), cs2(j,k-1,l), cs2(j,k,l), cs2(j,k+1,l), cs2(j,k+2,l), cs2R(j,k-1,l), cs2L(j,k,l))
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
! Do reconstruction using PPM interpolation, with the Mignone 2014 algorithm
! Reconstruct along the y-direction, and get states at the y-interfaces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PPM2014_recony
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
				CALL PPM_2014 (dy2(k-2), dy2(k-1), dy2(k), dy2(k+1), dy2(k+2), prim2(i,j,k-2,l), prim2(i,j,k-1,l), prim2(i,j,k,l), prim2(i,j,k+1,l), prim2(i,j,k+2,l), primR2(i,j,k-1,l), primL2(i,j,k,l))
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
				CALL PPM_2014 (dy2(k-2), dy2(k-1), dy2(k), dy2(k+1), dy2(k+2), cs2(j,k-2,l), cs2(j,k-1,l), cs2(j,k,l), cs2(j,k+1,l), cs2(j,k+2,l), cs2R(j,k-1,l), cs2L(j,k,l))
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
! Do reconstruction using PPM interpolation, with the original Colella algorithm
! Reconstruct along the z-direction, and get states at the z-interfaces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PPM1984_reconz
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
				CALL PPM_1984 (dz2(l-2), dz2(l-1), dz2(l), dz2(l+1), dz2(l+2), prim2(i,j,k,l-2), prim2(i,j,k,l-1), prim2(i,j,k,l), prim2(i,j,k,l+1), prim2(i,j,k,l+2), primR2(i,j,k,l-1), primL2(i,j,k,l))
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
				CALL PPM_1984 (dz2(l-2), dz2(l-1), dz2(l), dz2(l+1), dz2(l+2), cs2(j,k,l-2), cs2(j,k,l-1), cs2(j,k,l), cs2(j,k,l+1), cs2(j,k,l+2), cs2R(j,k,l-1), cs2L(j,k,l))
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do reconstruction using PPM interpolation, with the Mignone 2011 algorithm
! Reconstruct along the z-direction, and get states at the z-interfaces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PPM2011_reconz
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
				CALL PPM_2011 (dz2(l-2), dz2(l-1), dz2(l), dz2(l+1), dz2(l+2), prim2(i,j,k,l-2), prim2(i,j,k,l-1), prim2(i,j,k,l), prim2(i,j,k,l+1), prim2(i,j,k,l+2), primR2(i,j,k,l-1), primL2(i,j,k,l))
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
				CALL PPM_2011 (dz2(l-2), dz2(l-1), dz2(l), dz2(l+1), dz2(l+2), cs2(j,k,l-2), cs2(j,k,l-1), cs2(j,k,l), cs2(j,k,l+1), cs2(j,k,l+2), cs2R(j,k,l-1), cs2L(j,k,l))
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do reconstruction using PPM interpolation, with the Mignone 2014 algorithm
! Reconstruct along the z-direction, and get states at the z-interfaces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PPM2014_reconz
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
				CALL PPM_2014 (dz2(l-2), dz2(l-1), dz2(l), dz2(l+1), dz2(l+2), prim2(i,j,k,l-2), prim2(i,j,k,l-1), prim2(i,j,k,l), prim2(i,j,k,l+1), prim2(i,j,k,l+2), primR2(i,j,k,l-1), primL2(i,j,k,l))
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
				CALL PPM_2014 (dz2(l-2), dz2(l-1), dz2(l), dz2(l+1), dz2(l+2), cs2(j,k,l-2), cs2(j,k,l-1), cs2(j,k,l), cs2(j,k,l+1), cs2(j,k,l+2), cs2R(j,k,l-1), cs2L(j,k,l))
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate cell average values to interface values using PPM 
! Assumed non-uniform gridding. Using the original PPM algorithm
! See Colella 1984. No steepening or flattening is performed.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PPM_1984 (dm2, dm1, dc, dp1, dp2, vm2, vm1, vc, vp1, vp2, vm_out, vp_out)
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

! Get coefficient !
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Getting interface primitive variables using the PPM method
! with the variants proposed by Mignone 2011
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PPM_2011 (dm2, dm1, dc, dp1, dp2, vm2, vm1, vc, vp1, vp2, vm_out, vp_out)
USE DEFINITION
IMPLICIT NONE

! The input into the subroutine, including conservative variable and input flux function !
REAL*8, INTENT (IN) :: dm2, dm1, dc, dp1, dp2
REAL*8, INTENT (IN) :: vm2, vm1, vc, vp1, vp2

! The output of the subroutine, the flux at cell boundary !
REAL*8, INTENT (OUT) :: vp_out, vm_out

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Getting interface primitive variables using the PPM method
! with the variants proposed by Mignone 2014
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PPM_2014 (dm2, dm1, dc, dp1, dp2, vm2, vm1, vc, vp1, vp2, vm_out, vp_out)
USE DEFINITION
IMPLICIT NONE

! The input into the subroutine, including conservative variable and input flux function !
REAL*8, INTENT (IN) :: dm2, dm1, dc, dp1, dp2
REAL*8, INTENT (IN) :: vm2, vm1, vc, vp1, vp2

! The output of the subroutine, the flux at cell boundary !
REAL*8, INTENT (OUT) :: vp_out, vm_out

END SUBROUTINE

END MODULE
