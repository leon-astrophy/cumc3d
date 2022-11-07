!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine update the pressure profile and their deriative    !
! once the density profile is being updated through rungekutta time  !
! evolution. It is being used in every time step, do not confused it !
! with subroutine GETRHOEOSRTOP                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDPRESSURE
USE OMP_LIB
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k, l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We do the DM case first !
! This would be done only if users wants DM component !
IF (RUNDM_flag) THEN
	DO CONCURRENT (j = nx_min_1:nx_part_1, k = ny_min_1:ny_part_1, l = nz_min_1:nz_part_1) 	
		p1 (j,k,l) = kgas1 * prim1(irho1,j,k,l) ** ggas1
		dpdrho1 (j,k,l) = kgas1 * ggas1 * prim1(irho1,j,k,l) ** (ggas1 - 1.0D0)
		dpdeps1 (j,k,l) = 0.0D0
	END DO

	! Copy to boundary !
	CALL BOUNDARY1D_DM(p1, even, part)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! The following steps are more or less similar , so no repeat 
DO CONCURRENT (j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2) 
	prim2(itau2,j,k,l) = prim2(irho2,j,k,l) * epsilon2(j,k,l) * (ggas2 - 1.0E0_DP) 
	dpdrho2 (j,k,l) = epsilon2(j,k,l) * (ggas2 - 1.0E0_DP)
	dpdeps2 (j,k,l) = prim2(irho2,j,k,l) * (ggas2 - 1.0E0_DP)
END DO

! Copy to boundary !
CALL BOUNDARY1D_NM(prim2(itau2,:,:,:), even, part)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

	real(DP) function large_pressure(x)
	implicit none
	real(DP) :: x
	large_pressure = x*SQRT(x**2 + 1.0D0)*(2.0D0*x**2 - 3.0D0) + 3.0D0*log(x + SQRT(x**2 + 1.0D0))
	end function

	real(DP) function small_pressure(x)
	implicit none
	real(DP) :: x
	small_pressure = 1.6D0*x**5 - (4.0D0/7.0D0)*x**7 + (1.0D0/3.0D0)*x**9 - (5.0D0/2.2D1)*x**11 & 
			+ (3.5D1/2.08D2)*x**13 - (2.1D1/1.6D2)*x**15 + (2.31D2/2.176D3)*x**17
	end function

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the specific internal energy !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EOSPRESSURE (rho_in, eps_in, type_in, p_out)
USE DEFINITION
IMPLICIT NONE

! Input type !
INTEGER, INTENT (IN) :: type_in

! Input density !
REAL (DP), INTENT (IN) :: rho_in, eps_in

! Output value !
REAL (DP), INTENT (OUT) :: p_out

! For DM/NM, choose by type !
IF(type_in == 1) THEN
	p_out = kgas1 * rho_in ** (ggas1)
ELSEIF(type_in == 2) THEN
	p_out = rho_in * eps_in ** (ggas2 - 1.0D0)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

	real(DP) function large_pressure(x)
	implicit none
	real(DP) :: x
	large_pressure = x*SQRT(x**2 + 1.0D0)*(2.0D0*x**2 - 3.0D0) + 3.0D0*log(x + SQRT(x**2 + 1.0D0))
	end function

	real(DP) function small_pressure(x)
	implicit none
	real(DP) :: x
	small_pressure = 1.6D0*x**5 - (4.0D0/7.0D0)*x**7 + (1.0D0/3.0D0)*x**9 - (5.0D0/2.2D1)*x**11 & 
			+ (3.5D1/2.08D2)*x**13 - (2.1D1/1.6D2)*x**15 + (2.31D2/2.176D3)*x**17
	end function

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the pressure gradient !                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDGRADP
USE DEFINITION
USE RIEMANN_MODULE
IMPLICIT NONE

! Integer !
INTEGER :: j, k, l, p

! Loop !
DO CONCURRENT (j = nx_min_2:nx_part_2, k = ny_min_2:ny_part_2, l = nz_min_2:nz_part_2) 
	dp_2 (x_dir,j,k,l) = (primL2(itau2,x_dir,j,k,l) - primR2(itau2,x_dir,j-1,k,l))/dr2(j)
	IF(n_dim > 1) THEN
		dp_2 (y_dir,j,k,l) = (primL2(itau2,y_dir,j,k,l) - primR2(itau2,y_dir,j,k-1,l))/dy2
	END IF
	IF(n_dim > 2) THEN
		dp_2 (z_dir,j,k,l) = (primL2(itau2,z_dir,j,k,l) - primR2(itau2,z_dir,j,k,l-1))/dz2
	END IF
END DO

END SUBROUTINE