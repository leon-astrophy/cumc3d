!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine update the pressure profile and their deriative    !
! once the density profile is being updated through rungekutta time  !
! evolution. It is being used in every time step, do not confused it !
! with subroutine GETRHOEOSRTOP                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDPRESSURE
USE RIEMANN_MODULE
USE HELMEOS_MODULE
USE DEFINITION
use nuceos_module
use ecap_module
IMPLICIT NONE

! Integer parameter !
INTEGER :: j, k, flag_eostable

! Dummy variables !
REAL (DP) :: dummy, dxdrho

! dummies for nuc eos
real (DP) :: rho_in
real (DP) :: p_out
real (DP) :: dpdrho_out
real (DP) :: dpdeps_out, depsdt_out
real (DP) :: mu_e_loc, mu_p_loc, mu_n_loc, mu_nu_loc, eta_loc

! dummies for EOS call
real*8 eosdummy
integer keyerr,keytemp
! end dummies for EOS call

keytemp = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We do the DM case first !
! This would be done only if users wants DM component !
IF (DM_flag == 1) THEN
	IF (fermieosdm_flag == 1) Then
		DO k = 1, length_step_z_1
			DO j = 1, length_step_r_1
				CALL FERMIMO (dlfmmo1, rho1 (j,k), 1)
				CALL FINDDXDRHO (dxdrho, rho1 (j,k), 1)
				IF (dlfmmo1 <= 1.0E-2_DP) THEN
					p1 (j,k) = a_max1*small_pressure(dlfmmo1)
				ELSE
					p1 (j,k) = a_max1*large_pressure(dlfmmo1)
				END IF

				! Derivative !
				dpdrho1 (j,k) = a_max1*dxdrho*dpdx(dlfmmo1)
				dpdepsilon1 (j,k) = 0.0E0_DP
			END DO
		END DO
	ELSEIF (polyeosdm_flag == 1) THEN
		p1 (:,:) = k_1 * rho1(:,:) ** gamma1
		dpdrho1 (:,:) = k_1 * gamma1 * rho1(:,:) ** (gamma1 - 1.0D0)
		dpdepsilon1 (:,:) = 0.0D0
	END IF
	CALL BOUNDARY1D_DMFULL (p1, even)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! The following steps are more or less similar , so no repeat !
IF (fermieosnm_flag == 1) THEN
	DO k = 1, length_step_z_2
		DO j = 1, length_step_r_2
			CALL FERMIMO (dlfmmo2, rho2 (j,k), 2)
			CALL FINDDXDRHO (dxdrho, rho2 (j,k), 2)
			IF (dlfmmo2 <= 1.0E-2_DP) THEN
				p2 (j,k) = a_max2*small_pressure(dlfmmo2)
			ELSE
				p2 (j,k) = a_max2*large_pressure(dlfmmo2)
			END IF

			! Derivative !
			dpdrho2 (j,k) = a_max2*dxdrho*dpdx(dlfmmo2)
		END DO
	END DO
ELSEIF (polyeosnm_flag == 1) THEN	
	IF(nm_epsilon == 0) THEN
		p2 (:,:) = k_2 * rho2(:,:) ** gamma2
		dpdrho2 (:,:) = k_2 * gamma2 * rho2(:,:) ** (gamma2 - 1.0D0)
		dpdepsilon2 (:,:) = 0.0D0
	ELSE
		p2 (:,:) = rho2(:,:) * epsilon2(:,:) * (gamma2 - 1.0E0_DP) 
		dpdrho2 (:,:) = epsilon2(:,:) * (gamma2 - 1.0E0_DP)
		dpdepsilon2 (:,:) = rho2(:,:) * (gamma2 - 1.0E0_DP)
	END IF
ELSEIF (helmeos_flag == 1) THEN
	DO k = 1, length_step_z_2
		DO j = 1, length_step_r_2
			CALL HELMEOS_RtoP(rho2 (j,k), temp2(j,k), abar2(j,k), zbar2(j,k), ye2(j,k), p2 (j,k), dpdrho2 (j,k), dpdepsilon2 (j,k))
		END DO
	END DO
ELSEIF (nuceos_flag == 1) THEN
	DO k = 1, length_step_z_2
		DO j = 1, length_step_r_2
	   		if(rho2(j,k) < rho2_a) THEN
				rho2(j,k) = rho2_a	   
			END IF
	   		rho_in = rho2(j,k) * 6.171D17
	   		call nuc_eos_full(rho_in,temp2(j,k),ye2(j,k),eosdummy,&
                              p_out,eosdummy,eosdummy,depsdt_out,&
                              dpdeps_out,dpdrho_out,eosdummy,&
                              eosdummy,eosdummy,eosdummy,eosdummy,&
                              eosdummy,mu_e_loc,mu_n_loc,mu_p_loc,&
							  mu_nu_loc,keytemp,keyerr,eos_rf_prec)
	   		p2(j,k) = p_out * 1.8005D-39
	   		dpdrho2(j,k) = dpdrho_out / 9.0D20
	   		dpdepsilon2(j,k) = dpdeps_out / 6.171D17
	   		depsdt_loc(j,k) = depsdt_out / 9.0D20

	   		if(NuStress_flag == 1 .and. bounce_flag == 0) then
	      		if(rho2(j,k) >= 3.2410D-6) then
		 			mu_nu_loc = mu_e_loc + mu_p_loc - mu_n_loc 
		 			eta_loc = MAX(mu_nu_loc / temp2(j,k), 0.0D0)
		 			p2nu(j,k) = 1.50985D-15 * 4.0D0 * pi / 3.0D0 * temp2(j,k)**4 * Fermi3(eta_loc)
		 			p2(j,k) = p2(j,k) + p2nu(j,k)
                 	dpdepsilon2(j,k) = dpdepsilon2(j,k) + &
		     		p2nu(j,k) / temp2(j,k) * (4.0D0 - eta_loc * Fermi3der(eta_loc) / Fermi3(eta_loc)) / depsdt_loc(j,k)
	      		endif
	   		endif
		END DO
	END DO
END IF

! Copy to boundary !
CALL BOUNDARY1D_NMFULL (p2, even)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

	real(DP) function large_pressure(x)
	implicit none
	real(DP) :: x
	large_pressure = x*SQRT(x**2 + 1.0D0)*(2.0D0*x**2 - 3.0D0) + 3.0D0*log(x + SQRT(x**2 + 1.0D0))
	end function

	real(DP) function dpdx(x)
	implicit none
	real(DP) :: x
	dpdx = 8.0D0*x**4/SQRT(x**2 + 1.0D0)
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
SUBROUTINE FINDDPDR
USE DEFINITION
USE RIEMANN_MODULE
IMPLICIT NONE

! Integer !
INTEGER :: j, k

! Loop !
IF(lapse_flag == 1) THEN
	DO j = 1, length_step_r_part_2
		dp2dr (j,:) = (ap2L(j,:)*p2L(j,:) - ap2R(j-1,:)*p2R(j-1,:))/dr2(j)
	END DO
ELSE
	DO j = 1, length_step_r_part_2
		dp2dr (j,:) = (p2L(j,:) - p2R(j-1,:))/dr2(j)
	END DO
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the pressure gradient !                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDDPDZ
USE DEFINITION
USE RIEMANN_MODULE
IMPLICIT NONE

! Integer !
INTEGER :: j, k

! Loop !
IF(lapse_flag == 1) THEN
	DO k = length_step_z_min_part_2, length_step_z_part_2
		dp2dz (:,k) = (ap2L(:,k)*p2L(:,k) - ap2R(:,k-1)*p2R(:,k-1))/dx2_z
	END DO
ELSE 
	DO k = length_step_z_min_part_2, length_step_z_part_2
		dp2dz (:,k) = (p2L(:,k) - p2R(:,k-1))/dx2_z
	END DO
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the specific internal energy !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EOSPRESSURE (den, pre, eps, type)
USE DEFINITION
IMPLICIT NONE

! Input type !
INTEGER, INTENT (IN) :: type

! Input density !
REAL (DP), INTENT (IN) :: den, eps

! Output value !
REAL (DP), INTENT (OUT) :: pre

! Extra patches by Ivan !
REAL (DP) :: p_poly, eps_poly, eps_thermal

! Fermi-momentum !
REAL (DP) :: fermi

! For DM Output !
IF(type == 1) THEN
	IF(fermieosdm_flag == 1) THEN
		CALL FERMIMO (fermi, den, 1)
		IF (fermi<=1.0E-2_DP) THEN
			pre = a_max1*small_pressure(fermi)
		ELSE
			pre = a_max1*large_pressure(fermi)
		END IF
	ELSE
		pre = k_1 * den ** (gamma1 - 1.0E0_DP) / (gamma1 - 1.0E0_DP)
	END IF

! For NM !
ELSEIF(type == 2) THEN
	IF(fermieosnm_flag == 1) THEN
		CALL FERMIMO (fermi, den, 2)
		IF (fermi<=1.0E-2_DP) THEN
			pre = a_max2*small_pressure(fermi)
		ELSE
			pre = a_max2*large_pressure(fermi)
		END IF
	ELSE
		IF(nm_epsilon == 0) THEN
			pre = k_2 * den ** (gamma2 - 1.0E0_DP) / (gamma2 - 1.0E0_DP)
		ELSE
			pre = pre/den/(gamma2 - 1.0D0)
		END IF
	END IF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

	real(DP) function large_pressure(x)
	implicit none
	real(DP) :: x
	large_pressure = x*SQRT(x**2 + 1.0D0)*(2.0D0*x**2 - 3.0D0) + 3.0D0*log(x + SQRT(x**2 + 1.0D0))
	end function

	real(DP) function dpdx(x)
	implicit none
	real(DP) :: x
	dpdx = 8.0D0*x**4/SQRT(x**2 + 1.0D0)
	end function

	real(DP) function small_pressure(x)
	implicit none
	real(DP) :: x
	small_pressure = 1.6D0*x**5 - (4.0D0/7.0D0)*x**7 + (1.0D0/3.0D0)*x**9 - (5.0D0/2.2D1)*x**11 & 
			+ (3.5D1/2.08D2)*x**13 - (2.1D1/1.6D2)*x**15 + (2.31D2/2.176D3)*x**17
	end function

END SUBROUTINE