!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the soundspeed of DM or NM !                              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SOUNDSPEED
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k, l

! For DM !
IF (DM_flag) THEN
	DO CONCURRENT (j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1)
		cs1(j,k,l) = sqrt(dpdrho1(j,k,l)+dpdeps1(j,k,l)*p1(j,k,l)/prim1(irho1,j,k,l)**(2.0E0_DP))
	END DO
	CALL BOUNDARY1D_DM (cs1, even, part)
ENDIF

! For NM !
DO CONCURRENT (j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2)
	cs2(j,k,l) = sqrt(dpdrho2(j,k,l)+dpdeps2(j,k,l)*prim2(itau2,j,k,l)/prim2(irho2,j,k,l)**(2.0E0_DP))
END DO

! Boundary !
CALL BOUNDARY1D_NM (cs2, even, part)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the speed of sound !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EOSSOUNDSPEED(p_in, rho_in, eps_in, cs_out, type_in)
USE DEFINITION
IMPLICIT NONE

! Input type !
INTEGER, INTENT(IN) :: type_in

! Input density !
REAL (DP), INTENT (IN) :: p_in, rho_in, eps_in

! Output value !
REAL (DP), INTENT (OUT) :: cs_out

! Local variables !
REAL (DP) :: dpdden, dpdeps

! We do the DM case first !
IF (type_in == 1) Then
	dpdden = kgas1 * ggas1 * rho_in ** (ggas1 - 1.0D0)
	dpdeps = 0.0D0
	cs_out = sqrt(dpdden+dpdeps*p_in/rho_in**(2.0E0_DP))

! For NM !
ELSEIF(type_in == 2) THEN
	dpdden = eps_in * (ggas2 - 1.0E0_DP)
	dpdeps = rho_in * (ggas2 - 1.0E0_DP)
	cs_out = sqrt(dpdden+dpdeps*p_in/rho_in**(2.0E0_DP))
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

	real(DP) function dpdx(x)
	implicit none
	real(DP) :: x
	dpdx = 8.0D0*x**4/SQRT(x**2 + 1.0D0)
	end function

END SUBROUTINE