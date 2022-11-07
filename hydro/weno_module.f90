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
REAL (DP), PARAMETER :: smallpara = 1.0D-40

! The C constants
REAL (DP), DIMENSION (-1 : 2, 0 : 2) :: c

! The D and tilde-D constants
REAL (DP), DIMENSION (0 : 2) :: d, td

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do reconstruction using WENO interpolation, but along the vertical directions !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE WENO_Reconstruct
USE OMP_LIB
USE RIEMANN_MODULE
USE DEFINITION  
IMPLICIT NONE

! Integer !
INTEGER :: j, k, l, i, p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP PARALLEL SHARED(prim2, primL2, primR2, cs2, cs2L, cs2R, eps2L, eps2R, prim1, primL1, primR1, cs1, cs1L, cs1R, eps1L, eps1R, p1L, p1R)

! Do the same for DM !
IF(DM_flag) THEN
	!$OMP DO SIMD COLLAPSE(3) SCHEDULE (STATIC)
    DO l = nz_min_1 - 1, nz_part_1 + 1
		DO k = ny_min_1 - 1, ny_part_1 + 1
         	DO j = nx_min_1 - 1, nx_part_1 + 1

		 		! reconstruct !
				DO i = imin1, imax1
					CALL WENO_MN (j, x_dir, prim1(i,j-2,k,l), prim1(i,j-1,k,l), prim1(i,j,k,l), prim1(i,j+1,k,l), prim1(i,j+2,k,l), primR1(i,x_dir,j-1,k,l), primL1(i,x_dir,j,k,l), dm_f)
				END DO

		    	! Get the epsilon at boundary !
				CALL EOSPRESSURE (primR1(irho1,x_dir,j-1,k,l), p1R(x_dir,j-1,k,l), eps1R(x_dir,j-1,k,l), dm_f)
				CALL EOSPRESSURE (primL1(irho1,x_dir,j,k,l), p1L(x_dir,j,k,l), eps1L(x_dir,j,k,l), dm_f)
				CALL EOSEPSILON (primR1(irho1,x_dir,j-1,k,l), p1R(x_dir,j-1,k,l), eps1R(x_dir,j-1,k,l), dm_f)
				CALL EOSEPSILON (primL1(irho1,x_dir,j,k,l), p1L(x_dir,j,k,l), eps1L(x_dir,j,k,l), dm_f)
				CALL EOSSOUNDSPEED (p1R(x_dir,j-1,k,l), primR1(irho1,x_dir,j-1,k,l), eps1R(x_dir,j-1,k,l), cs1R(x_dir,j-1,k,l), dm_f)
				CALL EOSSOUNDSPEED (p1L(x_dir,j,k,l), primL1(irho1,x_dir,j,k,l), eps1L(x_dir,j,k,l), cs1L(x_dir,j,k,l), dm_f)

				! y-sweep !
				IF(n_dim > 1) THEN
					DO i = imin1, imax1
						CALL WENO_JS (k, prim1(i,j,k-2,l), prim1(i,j,k-1,l), prim1(i,j,k,l), prim1(i,j,k+1,l), prim1(i,j,k+2,l), primR1(i,y_dir,j,k-1,l), primL1(i,y_dir,j,k,l))
					END DO
					CALL EOSPRESSURE (primR1(irho1,y_dir,j,k-1,l), p1R(y_dir,j,k-1,l), eps1R(y_dir,j,k-1,l), dm_f)
					CALL EOSPRESSURE (primL1(irho1,y_dir,j,k,l), p1L(y_dir,j,k,l), eps1L(y_dir,j,k,l), dm_f)
					CALL EOSEPSILON (primR1(irho1,y_dir,j,k-1,l), p1R(y_dir,j,k-1,l), eps1R(y_dir,j,k-1,l), dm_f)
					CALL EOSEPSILON (primL1(irho1,y_dir,j,k,l), p1L(y_dir,j,k,l), eps1L(y_dir,j,k,l), dm_f)
					CALL EOSSOUNDSPEED (p1R(y_dir,j,k-1,l), primR1(irho1,y_dir,j,k-1,l), eps1R(y_dir,j,k-1,l), cs1R(y_dir,j,k-1,l), dm_f)
					CALL EOSSOUNDSPEED (p1L(y_dir,j,k,l), primL1(irho1,y_dir,j,k,l), eps1L(y_dir,j,k,l), cs1L(y_dir,j,k,l), dm_f)
				END IF

				! z-sweep !
				IF(n_dim > 1) THEN
					DO i = imin1, imax1
						CALL WENO_JS (l, prim1(i,j,k,l-2), prim1(i,j,k,l-1), prim1(i,j,k,l), prim1(i,j,k,l+1), prim1(i,j,k,l+2), primR1(i,z_dir,j,k,l-1), primL1(i,z_dir,j,k,l))
					END DO
					CALL EOSPRESSURE (primR1(irho1,z_dir,j,k,l-1), p1R(z_dir,j,k,l-1), eps1R(z_dir,j,k,l-1), dm_f)
					CALL EOSPRESSURE (primL1(irho1,z_dir,j,k,l), p1L(z_dir,j,k,l), eps1L(z_dir,j,k,l), dm_f)
					CALL EOSEPSILON (primR1(irho1,z_dir,j,k,l-1), p1R(z_dir,j,k,l-1), eps1R(z_dir,j,k,l-1), dm_f)
					CALL EOSEPSILON (primL1(irho1,z_dir,j,k,l), p1L(z_dir,j,k,l), eps1L(z_dir,j,k,l), dm_f)
					CALL EOSSOUNDSPEED (p1R(z_dir,j,k,l-1), primR1(irho1,z_dir,j,k,l-1), eps1R(z_dir,j,k,l-1), cs1R(z_dir,j,k,l-1), dm_f)
					CALL EOSSOUNDSPEED (p1L(z_dir,j,k,l), primL1(irho1,z_dir,j,k,l), eps1L(z_dir,j,k,l), cs1L(z_dir,j,k,l), dm_f)
				END IF

         	END DO
		END DO
	END DO
	!$OMP END DO
END IF

! We first interpolate density for NM !
!$OMP DO SIMD COLLAPSE(3) SCHEDULE (STATIC)
DO l = nz_min_2 - 1, nz_part_2 + 1
	DO k = ny_min_2 - 1, ny_part_2 + 1
		DO j = nx_min_2 - 1, nx_part_2 + 1
   		
			! Reconstruct for the x-direction !
			DO i = imin2, imax2
				CALL WENO_MN (j, x_dir, prim2(i,j-2,k,l), prim2(i,j-1,k,l), prim2(i,j,k,l), prim2(i,j+1,k,l), prim2(i,j+2,k,l), primR2(i,x_dir,j-1,k,l), primL2(i,x_dir,j,k,l), nm_f)
			END DO
			IF(dual_energy) THEN
				eps2R(x_dir,j-1,k,l) = primR2(ieps2,x_dir,j-1,k,l)/primR2(irho2,x_dir,j-1,k,l)
				eps2L(x_dir,j,k,l) = primL2(ieps2,x_dir,j,k,l)/primL2(irho2,x_dir,j,k,l)
				CALL WENO_MN (j, x_dir, cs2(j-2,k,l), cs2(j-1,k,l), cs2(j,k,l), cs2(j+1,k,l), cs2(j+2,k,l), cs2R(x_dir,j-1,k,l), cs2L(x_dir,j,k,l), nm_f)
			ELSE
			   	CALL EOSEPSILON (primR2(irho2,x_dir,j-1,k,l), primR2(itau2,x_dir,j-1,k,l), eps2R(x_dir,j-1,k,l), nm_f)
			   	CALL EOSEPSILON (primL2(irho2,x_dir,j,k,l), primL2(itau2,x_dir,j,k,l), eps2L(x_dir,j,k,l), nm_f)
			   	CALL EOSSOUNDSPEED (primR2(itau2,x_dir,j-1,k,l), primR2(irho2,x_dir,j-1,k,l), eps2R(x_dir,j-1,k,l), cs2R(x_dir,j-1,k,l), nm_f)
			   	CALL EOSSOUNDSPEED (primL2(itau2,x_dir,j,k,l), primL2(irho2,x_dir,j,k,l), eps2L(x_dir,j,k,l), cs2L(x_dir,j,k,l), nm_f)
			END IF

			! y-sweep !
			IF(n_dim > 1) THEN
				DO i = imin2, imax2
					CALL WENO_JS (k, prim2(i,j,k-2,l), prim2(i,j,k-1,l), prim2(i,j,k,l), prim2(i,j,k+1,l), prim2(i,j,k+2,l), primR2(i,y_dir,j,k-1,l), primL2(i,y_dir,j,k,l))
				END DO
				IF(dual_energy) THEN
					eps2R(y_dir,j,k-1,l) = primR2(ieps2,y_dir,j,k-1,l)/primR2(irho2,y_dir,j,k-1,l)
					eps2L(y_dir,j,k,l) = primL2(ieps2,y_dir,j,k,l)/primL2(irho2,y_dir,j,k,l)
					CALL WENO_JS (k, cs2(j,k-2,l), cs2(j,k-1,l), cs2(j,k,l), cs2(j,k+1,l), cs2(j,k+2,l), cs2R(y_dir,j,k-1,l), cs2L(y_dir,j,k,l))
				ELSE
			   		CALL EOSEPSILON (primR2(irho2,y_dir,j,k-1,l), primR2(itau2,y_dir,j,k-1,l), eps2R(y_dir,j,k-1,l), nm_f)
			   		CALL EOSEPSILON (primL2(irho2,y_dir,j,k,l), primL2(itau2,y_dir,j,k,l), eps2L(y_dir,j,k,l), nm_f)
			   		CALL EOSSOUNDSPEED (primR2(itau2,y_dir,j,k-1,l), primR2(irho2,y_dir,j,k-1,l), eps2R(y_dir,j,k-1,l), cs2R(y_dir,j,k-1,l), nm_f)
			   		CALL EOSSOUNDSPEED (primL2(itau2,y_dir,j,k,l), primL2(irho2,y_dir,j,k,l), eps2L(y_dir,j,k,l), cs2L(y_dir,j,k,l), nm_f)
				END IF
			END IF

			! z-sweep !
			IF(n_dim > 2) THEN
				DO i = imin2, imax2
					CALL WENO_JS (l, prim2(i,j,k,l-2), prim2(i,j,k,l-2), prim2(i,j,k,l), prim2(i,j,k,l+1), prim2(i,j,k,l+2), primR2(i,z_dir,j,k,l-1), primL2(i,z_dir,j,k,l))
				END DO
				IF(dual_energy) THEN
					eps2R(z_dir,j,k,l-1) = primR2(ieps2,z_dir,j,k,l-1)/primR2(irho2,z_dir,j,k,l-1)
					eps2L(z_dir,j,k,l) = primL2(ieps2,z_dir,j,k,l)/primL2(irho2,z_dir,j,k,l)
					CALL WENO_JS (l, cs2(j,k,l-2), cs2(j,k,l-1), cs2(j,k,l), cs2(j,k,l+1), cs2(j,k,l+2), cs2R(z_dir,j,k,l-1), cs2L(z_dir,j,k,l))
				ELSE
			   		CALL EOSEPSILON (primR2(irho2,z_dir,j,k,l-1), primR2(itau2,z_dir,j,k,l-1), eps2R(z_dir,j,k,l-1), nm_f)
			   		CALL EOSEPSILON (primL2(irho2,z_dir,j,k,l), primL2(itau2,z_dir,j,k,l), eps2L(z_dir,j,k,l), nm_f)
			   		CALL EOSSOUNDSPEED (primR2(itau2,z_dir,j,k,l-1), primR2(irho2,z_dir,j,k,l-1), eps2R(z_dir,j,k,l-1), cs2R(z_dir,j,k,l-1), nm_f)
			   		CALL EOSSOUNDSPEED (primL2(itau2,z_dir,j,k,l), primL2(irho2,z_dir,j,k,l), eps2L(z_dir,j,k,l), cs2L(z_dir,j,k,l), nm_f)
				END IF
			END IF

		END DO
	END DO
END DO
!$OMP END DO

!$OMP END PARALLEL
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the WENO scheme for reconstructing the numerical flux at both the !
! left and right hand side located at the boundary cell. In this version, I !
! provide different WENO scheme that differ by their smoothness indicator   !
! I also include WENO scheme that use combination of high and low order     !
! polynominal as building block. Nonetheless, a monotonicity preserving     !
! limter option is provided so to make the solution to be MPW               !
! For details, please refer to the textbook with ISBN 3-540-65893-9         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE WENO_MN (i_in, dir_in, vm2, vm1, vc, vp1, vp2, vm_out, vp_out, type_in)
USE DEFINITION
IMPLICIT NONE

! Input integer !
INTEGER, INTENT(IN) :: i_in, type_in, dir_in

! Input conservative variables !
REAL (DP), INTENT (IN) :: vm2, vm1, vc, vp1, vp2

! The reconstructed states at cell boundary !
REAL (DP), INTENT (OUT) :: vm_out, vp_out

! Integer !
INTEGER :: i, j, k

! Tempeorary parameter !
REAL (DP) :: tau5, temp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Patches from FEST-3D !
REAL (DP), DIMENSION (1 : 3) :: P_weno !< polynomial approximation
REAL (DP), DIMENSION (1 : 3) :: B_weno !< smoothness factor
REAL (DP), DIMENSION (1 : 3) :: w_weno !< wieght
REAL (DP), DIMENSION (1 : 3) :: g_weno !< linear wieght

! Real variables !
REAL (DP) :: U11_weno
REAL (DP) :: U00_weno
REAL (DP) :: U21_weno
REAL (DP) :: U10_weno
REAL (DP) :: U01_weno
REAL (DP) :: U12_weno
REAL (DP) :: alpha12_weno
REAL (DP) :: alpha01_weno
REAL (DP) :: alpha10_weno
REAL (DP) :: alpha21_weno

! Temporal arrays !
REAL (DP), DIMENSION (i_in - 2 : i_in + 2) :: v_weno
REAL (DP), DIMENSION (i_in - 2 : i_in + 2) :: vol_weno

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Linear weight !
g_weno (1) = 1.0D0/10.0D0
g_weno (2) = 6.0D0/10.0D0
g_weno (3) = 3.0D0/10.0D0

! Assign temporal arrays !
v_weno (i_in - 2) = vm2
v_weno (i_in - 1) = vm1
v_weno (i_in) = vc
v_weno (i_in + 1) = vp1
v_weno (i_in + 2) = vp2

! Grid size !
IF(type_in == 1) THEN
	IF(dir_in == 1) THEN
		vol_weno (i_in - 2) = dr1(i_in - 2)
		vol_weno (i_in - 1) = dr1(i_in - 1)
		vol_weno (i_in) = dr1(i_in)
		vol_weno (i_in + 1) = dr1(i_in + 1)
		vol_weno (i_in + 2) = dr1(i_in + 2)
	ELSEIF(dir_in == 2) THEN
		vol_weno (i_in - 2) = dy1
		vol_weno (i_in - 1) = dy1
		vol_weno (i_in) = dy1
		vol_weno (i_in + 1) = dy1
		vol_weno (i_in + 2) = dy1
	ELSEIF(dir_in == 3) THEN
		vol_weno (i_in - 2) = dz1
		vol_weno (i_in - 1) = dz1
		vol_weno (i_in) = dz1
		vol_weno (i_in + 1) = dz1
		vol_weno (i_in + 2) = dz1
	END IF
ELSEIF(type_in == 2) THEN
	IF(dir_in == 1) THEN
		vol_weno (i_in - 2) = dr2(i_in - 2)
		vol_weno (i_in - 1) = dr2(i_in - 1)
		vol_weno (i_in) = dr2(i_in)
		vol_weno (i_in + 1) = dr2(i_in + 1)
		vol_weno (i_in + 2) = dr2(i_in + 2)
	ELSEIF(dir_in == 2) THEN
		vol_weno (i_in - 2) = dy2
		vol_weno (i_in - 1) = dy2
		vol_weno (i_in) = dy2
		vol_weno (i_in + 1) = dy2
		vol_weno (i_in + 2) = dy2
	ELSEIF(dir_in == 3) THEN
		vol_weno (i_in - 2) = dz2
		vol_weno (i_in - 1) = dz2
		vol_weno (i_in) = dz2
		vol_weno (i_in + 1) = dz2
		vol_weno (i_in + 2) = dz2
	END IF
END IF

! Assign alpha !
alpha12_weno = vol_weno(i_in + 2)/(vol_weno(i_in + 1) + vol_weno(i_in + 2))
alpha01_weno = vol_weno(i_in + 1)/(vol_weno(i_in) + vol_weno(i_in + 1))
alpha10_weno = vol_weno(i_in)/(vol_weno(i_in - 1) + vol_weno(i_in))
alpha21_weno = vol_weno(i_in - 1)/(vol_weno(i_in - 2) + vol_weno(i_in - 1))

! Assign slopes !
U01_weno = (1.0D0 - alpha01_weno)*v_weno(i_in) + alpha01_weno*v_weno(i_in + 1)
U12_weno = (1.0D0 - alpha12_weno)*v_weno(i_in + 1) + alpha12_weno*v_weno(i_in + 2)
U10_weno = (1.0D0 - alpha10_weno)*v_weno(i_in - 1) + alpha10_weno*v_weno(i_in)
U21_weno = (1.0D0 - alpha21_weno)*v_weno(i_in - 2) + alpha21_weno*v_weno(i_in - 1)
U00_weno = v_weno(i_in - 1) + (1.0D0 - alpha21_weno)*(v_weno(i_in - 1) - v_weno(i_in - 2))
U11_weno = v_weno(i_in + 1) + alpha12_weno*(v_weno(i_in + 1) - v_weno(i_in + 2))

! Assign polynominals !
P_weno(1) = ( 6.0D0*v_weno(i_in) - 1.0D0*U10_weno - 2.0D0*U00_weno)/3.0D0
P_weno(2) = (- 1.0D0*U10_weno + 2.0D0*v_weno(i_in) + 2.0D0*U01_weno)/3.0D0
P_weno(3) = ( 2.0D0*U01_weno + 2.0D0*v_weno(i_in + 1) - 1.0D0*U12_weno)/3.0D0

! Assign smoothness indicator !
B_weno(1) = (13.0D0/12.0D0)*(2.0D0*U10_weno - 2.0*U00_weno)**2 + & 
	    (1.0D0/4.0D0)*(4.0D0*v_weno(i_in) - 2.0D0*U10_weno - 2.0D0*U00_weno)**2
B_weno(2) = (13.0D0/12.0D0)*(2.0D0*U10_weno - 4.0D0*v_weno(i_in) + 2.0D0*U01_weno)**2 + &
	    (1.0D0/4.0D0)*(2.0D0*U01_weno - 2.0D0*U10_weno)**2
B_weno(3) = (13.0D0/12.0D0)*(2.0D0*U01_weno - 4.0D0*v_weno(i_in + 1) + 2.0D0*U12_weno)**2 + & 
	    (1.0D0/4.0D0)*(8.0D0*v_weno(i_in + 1) - 6.0D0*U01_weno - 2.0D0*U12_weno)**2

! WENO-Z !
tau5 = abs(B_weno(1) - B_weno(3))

! Weight !
w_weno(:) = g_weno(:)*(1.0D0 + (tau5/(smallpara + B_weno(:)))**2)

! Output !
vp_out = SUM(w_weno*P_weno)/SUM(w_weno)

! Do the same for the right states !
P_weno(1) = (6.0D0*v_weno(i_in) - 1.0D0*U01_weno - 2.0D0*U11_weno)/3.0D0
P_weno(2) = (-1.0D0*U01_weno + 2.0D0*v_weno(i_in) + 2.0D0*U10_weno)/3.0D0
P_weno(3) = (2.0D0*U10_weno + 2.0D0*v_weno(i_in - 1) - 1.0D0*U21_weno)/3.0D0

B_weno(1) = (13.0D0/12.0D0)*(2.0D0*U01_weno - 2.0D0*U11_weno)**2 + & 
	    (1.0D0/4.0D0)*(4.0D0*v_weno(i_in) - 2.0D0*U01_weno - 2.0D0*U11_weno)**2
B_weno(2) = (13.0D0/12.0D0)*(2.0D0*U01_weno - 4.0D0*v_weno(i_in) + 2.0D0*U10_weno)**2 + &
	    (1.0D0/4.0D0)*(-2.0D0*u01_weno + 2.0D0*U10_weno)**2
B_weno(3) = (13.0D0/12.0D0)*(2.0D0*U10_weno - 4.0D0*v_weno(i_in - 1) + 2.0D0*U21_weno)**2 + & 
	    (1.0D0/4.0D0)*(-6.0D0*U10_weno + 8.0D0*v_weno(i_in - 1) - 2.0D0*U21_weno)**2

tau5 = abs(B_weno(1) - B_weno(3))
w_weno(:) = g_weno(:)*(1.0D0 + (tau5/(smallpara + B_weno(:)))**2)
vm_out = SUM(w_weno*P_weno)/SUM(w_weno)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine reads in the constant for WENO reconstuction
! assuming uniform grid anywhere along the row/column
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GetWenoConst
implicit none

! Dummy variable
integer :: r

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Const C

c (-1, 0) = 1.1E1_DP / 6.0E0_DP
c (-1, 1) = - 7.0E0_DP / 6.0E0_DP
c (-1, 2) = 1.0E0_DP / 3.0E0_DP

c (0, 0) = 1.0E0_DP / 3.0E0_DP
c (0, 1) = 5.0E0_DP / 6.0E0_DP
c (0, 2) = - 1.0E0_DP / 6.0E0_DP

c (1, 0) = - 1.0E0_DP / 6.0E0_DP
c (1, 1) = 5.0E0_DP / 6.0E0_DP
c (1, 2) = 1.0E0_DP / 3.0E0_DP

c (2, 0) = 1.0E0_DP / 3.0E0_DP
c (2, 1) = - 7.0E0_DP / 6.0E0_DP
c (2, 2) = 1.1E1_DP / 6.0E0_DP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Const D

d (0) = 3.0E0_DP / 1.0E1_DP
d (1) = 3.0E0_DP / 5.0E0_DP
d (2) = 1.0E0_DP / 1.0E1_DP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Const tilde_D

DO r = 0, 2
    td (r) = d (2 - r)
END DO

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the WENO scheme for reconstructing the numerical flux at both the !
! left and right hand side located at the boundary cell. In this version, I !
! provide different WENO scheme that differ by their smoothness indicator   !
! I also include WENO scheme that use combination of high and low order     !
! polynominal as building block. Nonetheless, a monotonicity preserving     !
! limter option is provided so to make the solution to be MPW               !
! For details, please refer to the textbook with ISBN 3-540-65893-9         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE WENO_JS (i, vm2, vm1, vc, vp1, vp2, vm_out, vp_out)
USE DEFINITION
IMPLICIT NONE

! Input integer !
INTEGER, INTENT(IN) :: i

! The input into the subroutine, including conservative variable and input flux function !
REAL (DP), INTENT (IN) :: vm2, vm1, vc, vp1, vp2

! The output of the subroutine, the flux at cell boundary !
REAL (DP), INTENT (OUT) :: vm_out, vp_out

! Temporal arrays !
REAL (DP), DIMENSION (0 : 2) :: vrhs, vlhs

! For assigning weights !
REAL (DP), DIMENSION (0 : 2) :: alpha, talpha, omega, tomega, beta

! Temporal arrays !
REAL (DP), DIMENSION (i - 2 : i + 2) :: v

! Integer !
INTEGER :: j, r, s

! Tempeorary parameter !
REAL (DP) :: tau, temp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign temporal arrays !
v(i - 2) = vm2
v(i - 1) = vm1
v(i) = vc
v(i + 1) = vp1
v(i + 2) = vp2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
! We calculate the value of u at each grid by the following loop !
! Do the right cell boundary !
DO r = 0, 2
	vrhs (r) = 0.0E0_DP
		
	! We calculate the value of u at right boundary !
	DO j = 0, 2
		vrhs (r) = vrhs (r) + c (r, j) * v (i - r + j)
	END DO
END DO

! Do the left cell boundary !
DO r = 0, 2
	vlhs (r) = 0.0E0_DP

	! Do the same for left boundary !
	DO j = 0, 2
		vlhs (r) = vlhs (r) + c (r - 1, j) * v (i - r + j)
	END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! These are essential parameter for further construction of u !
beta (0) = (1.3E1_DP / 1.2E1_DP) * (v (i) - 2 * v (i + 1) + v (i + 2)) ** 2 &
		+ (1.0E0_DP / 4.0E0_DP) * (3 * v (i) - 4 * v (i + 1) + v (i + 2)) ** 2
beta (1) = (1.3E1_DP / 1.2E1_DP) * (v (i - 1) - 2 * v (i) + v (i + 1)) ** 2 &
		+ (1.0E0_DP / 4.0E0_DP) * (v (i - 1) - v (i + 1)) ** 2
beta (2) = (1.3E1_DP / 1.2E1_DP) * (v (i - 2) - 2 * v (i - 1) + v (i)) ** 2 &
		+ (1.0E0_DP / 4.0E0_DP) * (v (i - 2) - 4 * v (i - 1) + 3 * v (i)) ** 2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assigning tau for the WENO-Z corrections !
tau = abs (beta(0) -beta(2))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do WENO-Z weight reconstructions !
DO r = 0, 2
	alpha (r) = d (r) * (1.0D0 + (tau/(beta(r) + smallpara))**2)
END DO

temp = 0.0E0_DP
	
! The denominator in finding omega, a coefficient for the last step of reconstruction  !
DO s = 0, 2
	temp = temp + alpha (s)
END DO
	
! Find the omega !
DO r = 0, 2
	omega (r) = alpha (r) / temp
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Find the alpha, omega for the value of u at left grid boundary... !
DO r = 0, 2
	talpha (r) = td (r) * (1.0D0 + (tau/(beta(r) + smallpara))**2)
END DO

temp = 0.0E0_DP

DO s = 0, 2
	temp = temp + talpha (s)
END DO

DO r = 0, 2
	tomega (r) = talpha (r) / temp
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

vp_out = 0.0E0_DP
	
! u at the left boundary !
DO r = 0, 2

	! Original WENO !
	vp_out = vp_out + omega (r) * vrhs (r)

END DO

vm_out = 0.0E0_DP
	
! u at the right boundary !
DO r = 0, 2

	! Original WENO !
	vm_out = vm_out + tomega (r) * vlhs (r)	

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

END MODULE