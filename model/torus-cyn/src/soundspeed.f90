!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the soundspeed of DM or NM !                              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SOUNDSPEED
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k, l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For NM !
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC)
DO l = nz_min_2 - 3, nz_part_2 + 3
	DO k = ny_min_2 - 3, ny_part_2 + 3
		DO j = nx_min_2 - 3, nx_part_2 + 3
			cs2(j,k,l) = sqrt(dpdrho2(j,k,l)+dpdeps2(j,k,l)*prim2(itau2,j,k,l)/prim2(irho2,j,k,l)**(2.0D0))
		END DO
	END DO
END DO
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the speed of sound !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EOSSOUNDSPEED(p_in, rho_in, eps_in, cs_out)
USE DEFINITION
IMPLICIT NONE

! Input density !
REAL*8, INTENT (IN) :: p_in, rho_in, eps_in

! Output value !
REAL*8, INTENT (OUT) :: cs_out

! Local variables !
REAL*8 :: dpdden, dpdeps

! We do the DM case first !
dpdden = eps_in * (ggas2 - 1.0D0)
dpdeps = rho_in * (ggas2 - 1.0D0)
cs_out = sqrt(dpdden+dpdeps*p_in/rho_in**(2.0D0))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

	REAL*8 function dpdx(x)
	implicit none
	REAL*8 :: x
	dpdx = 8.0D0*x**4/SQRT(x**2 + 1.0D0)
	end function

END SUBROUTINE