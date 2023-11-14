!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the specific internal energy !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EOSEPSILON_NM (rho_in, p_in, eps_out)
!$acc routine seq
USE DEFINITION
IMPLICIT NONE

! Input density !
REAL*8, INTENT (IN) :: rho_in, p_in

! Output value !
REAL*8, INTENT (OUT) :: eps_out

! For DM Output !
eps_out = p_in/rho_in/(ggas - 1.0D0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

	REAL*8 function large_energy(x)
	!$acc routine seq
	implicit none
	REAL*8 :: x
	large_energy = 3.0D0*x*DSQRT(x**2 + 1.0D0)*(1.0D0 + 2.0D0*x**2) - 3.0D0*log(x + DSQRT(x**2 + 1.0D0)) - 8.0D0*x**3
	end function

	REAL*8 function small_energy(x)
	!$acc routine seq
	implicit none
	REAL*8 :: x
	small_energy = 8.0D0*x**3 + (1.2D1/5.0D0)*x**5 - (3.0D0/7.0D0)*x**7 + (1.0D0/6.0D0)*x**9 - (1.5D1/1.76D2)*x**11 & 
							 + (2.1D1/4.16D2)*x**13 - (2.1D1/6.40D2)*x**15 + (9.9D1/4.352D3)*x**17 - 8.0D0*x**3
	end function

END SUBROUTINE