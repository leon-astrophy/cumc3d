!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Riemann problem tests for 1D hydrodynamics 
! 1. Toro test 1
! 2. Toro test 2
! 3. Toro test 3
! 4. Toro test 4
! 5. Toro test 5
! 6. Two interacting shocks
! 7. SHU-Osher problem
! 8. Sedov explosion
! 9. Dai & Woodward (1994)
! 10. Brio-Wu MHD shock tube 
! 11. Seven discontinutieis
! 12. Magnetosonic rarefractions
! 13. Two fast shocks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Riemann_1d
USE DEFINITION
IMPLICIT NONE
INCLUDE "param.h"

! Integer parameter !
integer :: j

! Real parameter !
REAL*8 :: mu_sedov, e_total

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Toro Rieamnn problem test 1, outflow boundaries both side !
if(test_model == 1) then

  ggas = 1.4D0
  total_time = 0.25D0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/1,1,1,1,1,1/)

	do j=1,nx
		if(x(j) <= 0.5D0) THEN
			prim(irho,j,:,:) = 1.0D0
			prim(itau,j,:,:) = 1.0D0
		else
			prim(irho,j,:,:) = 0.125D0
			prim(itau,j,:,:) = 0.1D0
    endif
		epsilon(j,:,:) = prim(itau,j,:,:) / prim(irho,j,:,:) / (ggas-1.0D0)
	enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Toro modified sod test, outflow boundaries both side !
elseif(test_model == 2) then

  ggas = 1.4D0
  total_time = 0.15D0
  output_profiletime = total_time/50.0D0
	
	!boundary_flag = (/1,1,1,1,1,1/)

	do j=1,nx
		if(x(j) <= 0.5D0) THEN
			prim(ivx,j,:,:) = -2.0D0
		else
			prim(ivx,j,:,:) = 2.0D0
		endif
		prim(irho,j,:,:) = 1.0D0
		prim(itau,j,:,:) = 0.4D0
		epsilon(j,:,:) = prim(itau,j,:,:) / prim(irho,j,:,:) / (ggas-1.0D0)
	enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Toro 3rd riemann problem test, outflow boundaries both side !
elseif(test_model == 3) then

  ggas = 1.4D0
  total_time = 0.012D0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/1,1,1,1,1,1/)

	do j=1,nx
		if(x(j) <= 0.5D0) THEN
			prim(itau,j,:,:) = 1000.0D0
		else
			prim(itau,j,:,:) = 0.01D0
		endif
		prim(irho,j,:,:) = 1.0D0
		prim(ivx,j,:,:) = 0.0D0
		epsilon(j,:,:) = prim(itau,j,:,:) / prim(irho,j,:,:) / (ggas-1.0D0)
	enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Toro riemann problem test 4, outflow boundaries both side !
elseif(test_model == 4) then

  ggas = 1.4D0
  total_time = 0.035D0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/1,1,1,1,1,1/)

	do j=1,nx
		if(x(j) <= 0.5D0) THEN
			prim(itau,j,:,:) = 0.01D0
		else
			prim(itau,j,:,:) = 100.0D0
		endif
		prim(irho,j,:,:) = 1.0D0
		prim(ivx,j,:,:) = 0.0D0
		epsilon(j,:,:) = prim(itau,j,:,:) / prim(irho,j,:,:) / (ggas-1.0D0)
	enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Toro riemann problem test 5, outflow boundaries both side !
elseif(test_model == 5) then

  ggas = 1.4D0
  total_time = 0.035D0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/1,1,1,1,1,1/)

	do j=1,nx
		if(x(j) <= 0.5D0) THEN
			prim(ivx,j,:,:) = 19.5975D0
			prim(itau,j,:,:) = 460.894D0
		else
			prim(ivx,j,:,:) = -6.19633D0
			prim(itau,j,:,:) = 46.0950D0
		endif
		prim(irho,j,:,:) = 5.99924D0
		epsilon(j,:,:) = prim(itau,j,:,:) / prim(irho,j,:,:) / (ggas-1.0D0)
	enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! two interacting shocks, reflecting boundaries both side !
elseif(test_model == 6) then

  ggas = 1.4D0
  total_time = 0.038D0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/2,2,1,1,1,1/)

	do j=1,nx
		if(x(j) <= 0.1D0) THEN
			prim(itau,j,:,:) = 1.0D3
		elseif(x(j) >= 0.9D0) THEN
			prim(itau,j,:,:) = 1.0D2
		else
			prim(itau,j,:,:) = 1.0D-2
		endif
		prim(irho,j,:,:) = 1.0D0
		prim(ivx,j,:,:) = 0.0D0
		epsilon(j,:,:) = prim(itau,j,:,:) / prim(irho,j,:,:) / (ggas-1.0D0)
	enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The SHU-Osher problem, outflow boundaries both side !
elseif(test_model == 7) then

  ggas = 1.4D0
  total_time = 0.178D0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/1,1,1,1,1,1/)

	do j=1,nx
		if(x(j) <= 0.1D0) THEN
			prim(irho,j,:,:) = 3.857143D0
			prim(ivx,j,:,:) = 2.629369D0
			prim(itau,j,:,:) = 10.3333D0
		else
			prim(irho,j,:,:) = 1.0D0 + 0.2D0*DSIN(5.0D1*x(j))
			prim(ivx,j,:,:) = 0.0D0
			prim(itau,j,:,:) = 1.0D0
		endif
		epsilon(j,:,:) = prim(itau,j,:,:) / prim(irho,j,:,:) / (ggas-1.0D0)
	enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The sedov explosion problem, one side is reflecting, another side is outgoing !
elseif(test_model == 8) then

  ggas = 1.4D0
  total_time = 0.47D0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/2,1,1,1,1,1/)

	mu_sedov = DBLE(coordinate_flag + 1)

	do j=1,nx
		if(j <= 4) then
			prim(itau,j,:,:) = (3.0D0 * (ggas - 1.0D0) * 1.0D0) & 
												 / ((1.0D0 + mu_sedov) * pi * x(4) ** mu_sedov)
		else
			prim(itau,j,:,:) = 1.0D-5
		endif
		prim(irho,j,:,:) = 1.0D0
		prim(ivx,j,:,:) = 0.0D0
		epsilon(j,:,:) = prim(itau,j,:,:) / prim(irho,j,:,:) / (ggas-1.0D0)
	enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Dai & Woodward (1994) shock tube test !
elseif(test_model == 9) then

  ggas = (5.0d0/3.0d0)
  total_time = 0.15D0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/1,1,1,1,1,1/)

	do j=0,nx
		if(x(j) <= 0.5D0) THEN
			prim(irho,j,:,:) = 1.0d0
			prim(itau,j,:,:) = 1.0d0
			prim(iby,j,:,:) = 1.0d0
		else
			prim(irho,j,:,:) = 0.2d0
			prim(itau,j,:,:) = 0.1d0
			prim(iby,j,:,:) = 0.0d0
		endif
		prim(ibx,j,:,:) = 1.0D0
		epsilon(j,:,:) = prim(itau,j,:,:) / prim(irho,j,:,:) / (ggas-1.0D0)
	enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Brio-Wu (or Balsara) MHD shock tube !
elseif(test_model == 10) then

  ggas = 2.0D0
  total_time = 0.1D0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/1,1,1,1,1,1/)

	do j=0,nx
		if(x(j) <= 0.5D0) THEN
			prim(irho,j,:,:) = 1.0D0
      prim(itau,j,:,:) = 1.0D0
			prim(iby,j,:,:) = 1.0D0
		else
			prim(irho,j,:,:) = 0.125D0
      prim(itau,j,:,:) = 0.1D0
			prim(iby,j,:,:) = -1.0D0
		endif
		prim(ibx,j,:,:) = 0.75D0
		epsilon(j,:,:) = prim(itau,j,:,:) / prim(irho,j,:,:) / (ggas-1.0D0)
	enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Seven discontinuities test !
elseif(test_model == 11) then

	ggas = (5.0d0/3.0d0)
  total_time = 0.2D0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/1,1,1,1,1,1/)

	do j=0,nx
		if(x(j) <= 0.5D0) THEN
			prim(irho,j,:,:) = 1.08d0
			prim(itau,j,:,:) = 0.95d0
			prim(ivx,j,:,:) = 1.2d0
			prim(ivy,j,:,:) = 0.01d0
			prim(ivz,j,:,:) = 0.5d0
			prim(iby,j,:,:) = 3.6d0/DSQRT(4.0d0*pi)
		else
			prim(irho,j,:,:) = 1.0d0
			prim(itau,j,:,:) = 1.0d0
			prim(ivx,j,:,:) = 0.0d0
			prim(ivy,j,:,:) = 0.0d0
			prim(ivz,j,:,:) = 0.0d0
			prim(iby,j,:,:) = 4.0d0/DSQRT(4.0d0*pi)
		endif
		prim(ibx,j,:,:) = 2.0d0/DSQRT(4.0d0*pi)
		prim(ibz,j,:,:) = 2.0d0/DSQRT(4.0d0*pi)
		epsilon(j,:,:) = prim(itau,j,:,:) / prim(irho,j,:,:) / (ggas-1.0D0)
	enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Magnetosonic rarefractions !
elseif(test_model == 12) then

	ggas = (5.0d0/3.0d0)
  total_time = 0.1D0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/1,1,1,1,1,1/)

	do j=0,nx
		if(x(j) <= 0.5D0) THEN
			prim(ivx,j,:,:) = -1.0d0
		else
			prim(ivx,j,:,:) = 1.0d0
		endif
		prim(irho,j,:,:) = 1.0d0
		prim(itau,j,:,:) = 1.0d0
		prim(iby,j,:,:) = 1.0d0
		epsilon(j,:,:) = prim(itau,j,:,:) / prim(irho,j,:,:) / (ggas-1.0D0)
	enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Two fast shock !
elseif(test_model == 13) then

	ggas = (5.0d0/3.0d0)
  total_time = 0.03D0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/1,1,1,1,1,1/)

	do j=0,nx
		if(x(j) <= 0.5D0) THEN
			prim(ivx,j,:,:) = 36.87d0
			prim(ivy,j,:,:) = -0.155d0
			prim(ivz,j,:,:) = -0.0386d0
		else
			prim(ivx,j,:,:) = -36.87d0
			prim(ivy,j,:,:) = 0.0d0
			prim(ivz,j,:,:) = 0.0d0
		endif
		prim(irho,j,:,:) = 1.0d0
		prim(itau,j,:,:) = 1.0d0
		prim(ibx,j,:,:) = 4.0d0/DSQRT(4.0d0*pi)
		prim(iby,j,:,:) = 4.0d0/DSQRT(4.0d0*pi)
		prim(ibz,j,:,:) = 1.0d0/DSQRT(4.0d0*pi)
		epsilon(j,:,:) = prim(itau,j,:,:) / prim(irho,j,:,:) / (ggas-1.0D0)
	enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! error message !
	
else
	stop "no such test model"
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

endsubroutine
