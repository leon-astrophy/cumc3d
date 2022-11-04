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

! The C constants
REAL (DP), DIMENSION (-1 : 2, 0 : 2) :: c

! The D and tilde-D constants
REAL (DP), DIMENSION (0 : 2) :: d, td

contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine reads in the constant for WENO reconstuction
   ! assuming uniform grid anywhere along the row/column
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine GetConst
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

   end subroutine getconst

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the WENO scheme for reconstructing the numerical flux at both the !
! left and right hand side located at the boundary cell. In this version, I !
! provide different WENO scheme that differ by their smoothness indicator   !
! I also include WENO scheme that use combination of high and low order     !
! polynominal as building block. Nonetheless, a monotonicity preserving     !
! limter option is provided so to make the solution to be MPW               !
! For details, please refer to the textbook with ISBN 3-540-65893-9         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE WENO (i, vm2, vm1, vc, vp1, vp2, vm_out, vp_out)
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
	alpha (r) = d (r) * (1.0D0 + (tau/(beta(r) + smallpara)))
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