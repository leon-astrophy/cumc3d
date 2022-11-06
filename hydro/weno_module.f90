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
USE RIEMANN_MODULE
USE DEFINITION  
IMPLICIT NONE

! Integer !
INTEGER :: j, k, l, i, p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP PARALLEL SHARED(prim2, primL2, primR2, cs2, cs2L, cs2R, eps2L, eps2R, vf2xR, vf2xL, vf2yR, vf2yL, vf2zR, vf2zL, xF2, yF2, zF2, &
!$OMP prim1, primL1, primR1, cs1, cs1L, cs1R, eps1L, eps1R, vf1xR, vf1xL, vf1yR, vf1yL, vf1zR, vf1zL, xF1, yF1, zF1)

! Do the same for DM !
IF(RUNDM_flag) THEN
	!$OMP DO SIMD COLLAPSE(3) SCHEDULE (STATIC)
	DO k = ny_min_1, ny_part_1
      DO l = nz_min_1, nz_part_1
         DO j = nx_min_1 - 1, nx_part_1 + 1
			   	DO i = imin1, imax1
				   	CALL WENO_MN (j, x_dir, prim1(j-2,k,l,i), prim1(j-1,k,l,i), prim1(j,k,l,i), prim1(j+1,k,l,i), prim1(j+2,k,l,i), primR1(j-1,k,l,i,x_dir), primL1(j,k,l,i,x_dir), dm_f)
			   	END DO

		       	! Get the epsilon at boundary !
			   	CALL EOSPRESSURE (primR1(j-1,k,l,irho1,x_dir), p1R(j-1,k,l,x_dir), eps1R(j-1,k,l,x_dir), dm_f)
			   	CALL EOSPRESSURE (primL1(j,k,l,irho1,x_dir), p1L(j,k,l,x_dir), eps1L(j,k,l,x_dir), dm_f)
			   	CALL EOSEPSILON (primR1(j-1,k,l,irho1,x_dir), p1R(j-1,k,l,x_dir), eps1R(j-1,k,l,x_dir), dm_f)
			   	CALL EOSEPSILON (primL1(j,k,l,irho1,x_dir), p1L(j,k,l,x_dir), eps1L(j,k,l,x_dir), dm_f)
			   	CALL EOSSOUNDSPEED (p1R(j-1,k,l,x_dir), primR1(j-1,k,l,irho1,x_dir), eps1R(j-1,k,l,x_dir), cs1R(j-1,k,l,x_dir), dm_f)
			   	CALL EOSSOUNDSPEED (p1L(j,k,l,x_dir), primL1(j,k,l,irho1,x_dir), eps1L(j,k,l,x_dir), cs1L(j,k,l,x_dir), dm_f)

				! No reconstruction for frame velocity since we assume analytic continous form !
				If(movinggridnm_flag) THEN
					vf1xR(j) = vel1_max*xF1(j)/radius1
					vf1xL(j) = vel1_max*xF1(j)/radius1
				END IF

         	END DO
		END DO
	END DO
	!$OMP END DO
END IF

! We first interpolate density for NM !
!$OMP DO SIMD COLLAPSE(3) SCHEDULE (STATIC)
DO j = nx_min_2 - 1, nx_part_2 + 1
	DO k = ny_min_2 - 1, ny_part_2 + 1
   		DO l = nz_min_2 - 1, nz_part_2 + 1

			! Reconstruct for the x-direction !
			DO i = imin2, imax2
				CALL WENO_MN (j, x_dir, prim2(j-2,k,l,i), prim2(j-1,k,l,i), prim2(j,k,l,i), prim2(j+1,k,l,i), prim2(j+2,k,l,i), primR2(j-1,k,l,i,x_dir), primL2(j,k,l,i,x_dir), nm_f)
			END DO
			IF(dual_energy) THEN
				eps2R(j-1,k,l,x_dir) = primR2(j-1,k,l,ieps2,x_dir)/primR2(j-1,k,l,irho2,x_dir)
				eps2L(j,k,l,x_dir) = primL2(j,k,l,ieps2,x_dir)/primL2(j,k,l,irho2,x_dir)
				CALL WENO_MN (j, x_dir, cs2(j-2,k,l), cs2(j-1,k,l), cs2(j,k,l), cs2(j+1,k,l), cs2(j+2,k,l), cs2R(j-1,k,l,x_dir), cs2L(j,k,l,x_dir), nm_f)
			ELSE
			   	CALL EOSEPSILON (primR2(j-1,k,l,irho2,x_dir), primR2(j-1,k,l,itau2,x_dir), eps2R(j-1,k,l,x_dir), nm_f)
			   	CALL EOSEPSILON (primL2(j,k,l,irho2,x_dir), primL2(j,k,l,itau2,x_dir), eps2L(j,k,l,x_dir), nm_f)
			   	CALL EOSSOUNDSPEED (primR2(j-1,k,l,itau2,x_dir), primR2(j-1,k,l,irho2,x_dir), eps2R(j-1,k,l,x_dir), cs2R(j-1,k,l,x_dir), nm_f)
			   	CALL EOSSOUNDSPEED (primL2(j,k,l,itau2,x_dir), primL2(j,k,l,irho2,x_dir), eps2L(j,k,l,x_dir), cs2L(j,k,l,x_dir), nm_f)
			END IF

			! No reconstruction for frame velocity since we assume analytic continous form !
			If(movinggridnm_flag) THEN
				vf2xR(j) = vel2_max*xF2(j)/radius2
				vf2xL(j) = vel2_max*xF2(j)/radius2
			END IF

			! y-sweep !
			IF(n_dim > 1) THEN
				DO i = imin2, imax2
					CALL WENO_MN (k, y_dir, prim2(j,k-2,l,i), prim2(j,k-1,l,i), prim2(j,k,l,i), prim2(j,k+1,l,i), prim2(j,k+2,l,i), primR2(j,k-1,l,i,y_dir), primL2(j,k,l,i,y_dir), nm_f)
				END DO
				IF(dual_energy) THEN
					eps2R(j,k-1,l,y_dir) = primR2(j,k-1,l,ieps2,y_dir)/primR2(j,k-1,l,irho2,y_dir)
					eps2L(j,k,l,y_dir) = primL2(j,k,l,ieps2,y_dir)/primL2(j,k,l,irho2,y_dir)
					CALL WENO_MN (k, y_dir, cs2(j,k-2,l), cs2(j,k-1,l), cs2(j,k,l), cs2(j,k+1,l), cs2(j,k+2,l), cs2R(j,k-1,l,y_dir), cs2L(j,k,l,y_dir), nm_f)
				ELSE
			   		CALL EOSEPSILON (primR2(j,k-1,l,irho2,y_dir), primR2(j,k-1,l,itau2,y_dir), eps2R(j,k-1,l,y_dir), nm_f)
			   		CALL EOSEPSILON (primL2(j,k,l,irho2,y_dir), primL2(j,k,l,itau2,y_dir), eps2L(j,k,l,y_dir), nm_f)
			   		CALL EOSSOUNDSPEED (primR2(j,k-1,l,itau2,y_dir), primR2(j,k-1,l,irho2,y_dir), eps2R(j,k-1,l,y_dir), cs2R(j,k-1,l,y_dir), nm_f)
			   		CALL EOSSOUNDSPEED (primL2(j,k,l,itau2,y_dir), primL2(j,k,l,irho2,y_dir), eps2L(j,k,l,y_dir), cs2L(j,k,l,y_dir), nm_f)
				END IF

				! No reconstruction for frame velocity since we assume analytic continous form !
				If(movinggridnm_flag) THEN
					vf2yR(k) = vel2_max*yF2(k)/radius2
					vf2yL(k) = vel2_max*yF2(k)/radius2
				END IF
			END IF

			! z-sweep !
			IF(n_dim > 2) THEN
				DO i = imin2, imax2
					CALL WENO_MN (k, z_dir, prim2(j,k,l-2,i), prim2(j,k,l-1,i), prim2(j,k,l,i), prim2(j,k,l+1,i), prim2(j,k,l+2,i), primR2(j,k,l-1,i,z_dir), primL2(j,k,l,i,z_dir), nm_f)
				END DO
				IF(dual_energy) THEN
					eps2R(j,k,l-1,z_dir) = primR2(j,k,l-1,ieps2,z_dir)/primR2(j,k,l-1,irho2,z_dir)
					eps2L(j,k,l,z_dir) = primL2(j,k,l,ieps2,z_dir)/primL2(j,k,l,irho2,z_dir)
					CALL WENO_MN (k, z_dir, cs2(j,k,l-2), cs2(j,k,l-1), cs2(j,k,l), cs2(j,k,l+1), cs2(j,k,l+2), cs2R(j,k,l-1,z_dir), cs2L(j,k,l,z_dir), nm_f)
				ELSE
			   		CALL EOSEPSILON (primR2(j,k,l-1,irho2,z_dir), primR2(j,k,l-1,itau2,z_dir), eps2R(j,k,l-1,z_dir), nm_f)
			   		CALL EOSEPSILON (primL2(j,k,l,irho2,z_dir), primL2(j,k,l,itau2,z_dir), eps2L(j,k,l,z_dir), nm_f)
			   		CALL EOSSOUNDSPEED (primR2(j,k,l-1,itau2,z_dir), primR2(j,k,l-1,irho2,z_dir), eps2R(j,k,l-1,z_dir), cs2R(j,k,l-1,z_dir), nm_f)
			   		CALL EOSSOUNDSPEED (primL2(j,k,l,itau2,z_dir), primL2(j,k,l,irho2,z_dir), eps2L(j,k,l,z_dir), cs2L(j,k,l,z_dir), nm_f)
				END IF

				! No reconstruction for frame velocity since we assume analytic continous form !
				If(movinggridnm_flag) THEN
					vf2zR(l) = vel2_max*zF2(l)/radius2
					vf2zL(l) = vel2_max*zF2(l)/radius2
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

END MODULE