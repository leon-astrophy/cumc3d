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
! 11. Falle 1998 test AW
! 12. Falle 1998 test FS
! 13. Falle 1998 test SS
! 14. Falle 1998 test FR
! 15. Falle 1998 test SR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Riemann_1d
USE DEFINITION
IMPLICIT NONE
INCLUDE "param.h"

! Integer parameter !
integer :: j

! Real parameter !
REAL*8 :: mu_sedov

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Toro Rieamnn problem test 1, outflow boundaries both side !
if(test_model == 1) then

  ggas2 = 1.4D0
  total_time = 0.25D0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/1,1,1,1,1,1/)

	do j=1,nx_2
		if(x2(j) <= 0.5D0) THEN
			prim2(irho2,j,:,:) = 1.0D0
			prim2(itau2,j,:,:) = 1.0D0
		else
			prim2(irho2,j,:,:) = 0.125D0
			prim2(itau2,j,:,:) = 0.1D0
    endif
		epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2-1.0D0)
	enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Toro modified sod test, outflow boundaries both side !
elseif(test_model == 2) then

  ggas2 = 1.4D0
  total_time = 0.15D0
  output_profiletime = total_time/50.0D0
	
	!boundary_flag = (/1,1,1,1,1,1/)

	do j=1,nx_2
		if(x2(j) <= 0.5D0) THEN
			prim2(ivel2_x,j,:,:) = -2.0D0
		else
			prim2(ivel2_x,j,:,:) = 2.0D0
		endif
		prim2(irho2,j,:,:) = 1.0D0
		prim2(itau2,j,:,:) = 0.4D0
		epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2-1.0D0)
	enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Toro 3rd riemann problem test, outflow boundaries both side !
elseif(test_model == 3) then

  ggas2 = 1.4D0
  total_time = 0.012D0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/1,1,1,1,1,1/)

	do j=1,nx_2
		if(x2(j) <= 0.5D0) THEN
			prim2(itau2,j,:,:) = 1000.0D0
		else
			prim2(itau2,j,:,:) = 0.01D0
		endif
		prim2(irho2,j,:,:) = 1.0D0
		prim2(ivel2_x,j,:,:) = 0.0D0
		epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2-1.0D0)
	enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Toro riemann problem test 4, outflow boundaries both side !
elseif(test_model == 4) then

  ggas2 = 1.4D0
  total_time = 0.035D0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/1,1,1,1,1,1/)

	do j=1,nx_2
		if(x2(j) <= 0.5D0) THEN
			prim2(itau2,j,:,:) = 0.01D0
		else
			prim2(itau2,j,:,:) = 100.0D0
		endif
		prim2(irho2,j,:,:) = 1.0D0
		prim2(ivel2_x,j,:,:) = 0.0D0
		epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2-1.0D0)
	enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Toro riemann problem test 5, outflow boundaries both side !
elseif(test_model == 5) then

  ggas2 = 1.4D0
  total_time = 0.035D0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/1,1,1,1,1,1/)

	do j=1,nx_2
		if(x2(j) <= 0.5D0) THEN
			prim2(ivel2_x,j,:,:) = 19.5975D0
			prim2(itau2,j,:,:) = 460.894D0
		else
			prim2(ivel2_x,j,:,:) = -6.19633D0
			prim2(itau2,j,:,:) = 46.0950D0
		endif
		prim2(irho2,j,:,:) = 5.99924D0
		epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2-1.0D0)
	enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! two interacting shocks, reflecting boundaries both side !
elseif(test_model == 6) then

  ggas2 = 1.4D0
  total_time = 0.038D0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/2,2,1,1,1,1/)

	do j=1,nx_2
		if(x2(j) <= 0.1D0) THEN
			prim2(itau2,j,:,:) = 1.0D3
		elseif(x2(j) >= 0.9D0) THEN
			prim2(itau2,j,:,:) = 1.0D2
		else
			prim2(itau2,j,:,:) = 1.0D-2
		endif
		prim2(irho2,j,:,:) = 1.0D0
		prim2(ivel2_x,j,:,:) = 0.0D0
		epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2-1.0D0)
	enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The SHU-Osher problem, outflow boundaries both side !
elseif(test_model == 7) then

  ggas2 = 1.4D0
  total_time = 0.178D0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/1,1,1,1,1,1/)

	do j=1,nx_2
		if(x2(j) <= 0.1D0) THEN
			prim2(irho2,j,:,:) = 3.857143D0
			prim2(ivel2_x,j,:,:) = 2.629369D0
			prim2(itau2,j,:,:) = 10.3333D0
		else
			prim2(irho2,j,:,:) = 1.0D0 + 0.2D0*DSIN(5.0D1*x2(j))
			prim2(ivel2_x,j,:,:) = 0.0D0
			prim2(itau2,j,:,:) = 1.0D0
		endif
		epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2-1.0D0)
	enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The sedov explosion problem, one side is reflecting, another side is outgoing !
elseif(test_model == 8) then

  ggas2 = 1.4D0
  total_time = 0.47D0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/2,1,1,1,1,1/)

	mu_sedov = DBLE(coordinate_flag + 1)

	do j=1,nx_2
		if(j <= 4) then
			prim2(itau2,j,:,:) = (3.0D0 * (ggas2 - 1.0D0) * 1.0D0) & 
												 / ((1.0D0 + mu_sedov) * pi * x2(4) ** mu_sedov)
		else
			prim2(itau2,j,:,:) = 1.0D-5
		endif
		prim2(irho2,j,:,:) = 1.0D0
		prim2(ivel2_x,j,:,:) = 0.0D0
		epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2-1.0D0)
	enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Dai & Woodward (1994) !
elseif(test_model == 9) then

  ggas2 = (5.d0 / 3.d0)
  total_time = 0.2D0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/1,1,1,1,1,1/)

	do j=0,nx_2
		if(x2(j) <= 0.5D0) THEN
			prim2(irho2,j,:,:) = 1.08D0
      prim2(itau2,j,:,:) = 0.95D0
			prim2(ivel2_x,j,:,:) = 1.2d0
			prim2(ivel2_y,j,:,:) = 0.01d0
			prim2(ivel2_z,j,:,:) = 0.5d0
		else
			prim2(irho2,j,:,:) = 1.0D0
      prim2(itau2,j,:,:) = 1.0D0
			prim2(ivel2_x,j,:,:) = 0.0d0
			prim2(ivel2_y,j,:,:) = 0.0d0
			prim2(ivel2_z,j,:,:) = 0.0d0
		endif
		prim2(ibx,j,:,:) = 2.d0 / DSQRT(16*DATAN(1.d0))
		prim2(iby,j,:,:) = 4.0d0 / DSQRT(16*DATAN(1.d0))
		prim2(ibz,j,:,:) = 2.d0 / DSQRT(16*DATAN(1.d0))
		epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2-1.0D0)
	enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Brio-Wu MHD shock tube !
elseif(test_model == 10) then

  ggas2 = 2.0D0
  total_time = 0.12D0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/1,1,1,1,1,1/)

	do j=0,nx_2
		if(x2(j) <= 0.5D0) THEN
			prim2(irho2,j,:,:) = 1.0D0
      prim2(itau2,j,:,:) = 1.0D0
			prim2(iby,j,:,:) = 1.0D0
		else
			prim2(irho2,j,:,:) = 0.125D0
      prim2(itau2,j,:,:) = 0.1D0
			prim2(iby,j,:,:) = -1.0D0
		endif
		prim2(ibx,j,:,:) = 0.75D0
		epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2-1.0D0)
	enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Falle AW test !
elseif(test_model == 11) then

	ggas2 = (5.0d0/3.0d0)
  total_time = 5.0D0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/1,1,1,1,1,1/)

	do j=0,nx_2
		prim2(irho2,j,:,:) = 1.0D0
		prim2(itau2,j,:,:) = 1.0D0
		prim2(ivel2_x,j,:,:) = 0.0D0
		prim2(ivel2_y,j,:,:) = 1.0D0
		prim2(ivel2_z,j,:,:) = 1.0D0
		prim2(ibx,j,:,:) = 1.0d0
		prim2(iby,j,:,:) = 1.0D0
		prim2(ibz,j,:,:) = 0.0D0
		epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2-1.0D0)
	enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Falle FS test !
elseif(test_model == 12) then

	ggas2 = (5.0d0/3.0d0)
  total_time = 0.4D0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/1,1,1,1,1,1/)

	do j=0,nx_2
		if(x2(j) <= 0.3D0) THEN
			prim2(irho2,j,:,:) = 3.0D0
			prim2(itau2,j,:,:) = 16.33d0
			prim2(ivel2_x,j,:,:) = -0.732D0
      prim2(ivel2_y,j,:,:) = -1.333D0
			prim2(iby,j,:,:) = 2.309D0
		else
			prim2(irho2,j,:,:) = 1.0D0
			prim2(itau2,j,:,:) = 1.0d0
			prim2(ivel2_x,j,:,:) = -4.196D0
      prim2(ivel2_y,j,:,:) = 0.0D0
			prim2(iby,j,:,:) = 0.0D0
		endif
		prim2(ivel2_z,j,:,:) = 0.0D0
		prim2(ibx,j,:,:) = 3.0d0
		prim2(ibz,j,:,:) = 0.0D0
		epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2-1.0D0)
	enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Falle SS test !
elseif(test_model == 13) then

	ggas2 = (5.0d0/3.0d0)
  total_time = 0.5D0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/1,1,1,1,1,1/)

	do j=0,nx_2
		if(x2(j) <= 0.3D0) THEN
			prim2(irho2,j,:,:) = 1.368D0
			prim2(itau2,j,:,:) = 1.769d0
			prim2(ivel2_x,j,:,:) = 0.269D0
      prim2(ivel2_y,j,:,:) = 1.0D0
			prim2(iby,j,:,:) = 0.0D0
		else
			prim2(irho2,j,:,:) = 1.0D0
			prim2(itau2,j,:,:) = 1.0D0
			prim2(ivel2_x,j,:,:) = 0.0D0
      prim2(ivel2_y,j,:,:) = 0.0D0
			prim2(iby,j,:,:) = 1.0D0
		endif
		prim2(ivel2_z,j,:,:) = 0.0D0
		prim2(ibx,j,:,:) = 1.0d0			
		prim2(ibz,j,:,:) = 0.0D0
		epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2-1.0D0)
	enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Falle FR test !
elseif(test_model == 14) then

	ggas2 = (5.0d0/3.0d0)
  total_time = 0.1D0
  output_profiletime = total_time/50.0D0

	!boundary_flag = (/1,1,1,1,1,1/)

	do j=0,nx_2
		if(x2(j) <= 0.5D0) THEN
			prim2(irho2,j,:,:) = 1.0d0
			prim2(itau2,j,:,:) = 2.0d0
			prim2(ivel2_x,j,:,:) = 0.0D0
      prim2(ivel2_y,j,:,:) = 0.0D0
			prim2(iby,j,:,:) = 3.0D0
		else
			prim2(irho2,j,:,:) = 0.2641d0
			prim2(itau2,j,:,:) = 0.2175d0
			prim2(ivel2_x,j,:,:) = 3.6d0
      prim2(ivel2_y,j,:,:) = -2.551d0
			prim2(iby,j,:,:) = 0.0D0
		endif
		prim2(ivel2_z,j,:,:) = 0.0D0
		prim2(ibx,j,:,:) = 1.0d0			
		prim2(ibz,j,:,:) = 0.0D0
		epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2-1.0D0)
	enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Falle SR test !
elseif(test_model == 15) then

	ggas2 = (5.0d0/3.0d0)
  total_time = 0.3D0
  output_profiletime = total_time/50.0D0
	
	!boundary_flag = (/1,1,1,1,1,1/)

	do j=0,nx_2
		if(x2(j) <= 0.5D0) THEN
			prim2(irho2,j,:,:) = 1.0d0
			prim2(itau2,j,:,:) = 2.0d0
			prim2(ivel2_x,j,:,:) = 0.0D0
      prim2(ivel2_y,j,:,:) = 0.0D0
			prim2(iby,j,:,:) = 0.0D0
		else
			prim2(irho2,j,:,:) = 0.2d0
			prim2(itau2,j,:,:) = 0.1368d0
			prim2(ivel2_x,j,:,:) = 1.186d0
      prim2(ivel2_y,j,:,:) = 2.967d0
			prim2(iby,j,:,:) = 1.6405d0
		endif
		prim2(ivel2_z,j,:,:) = 0.0D0
		prim2(ibx,j,:,:) = 1.0d0
		prim2(ibz,j,:,:) = 0.0D0
		epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2-1.0D0)
	enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! error message !
else
	stop "no such test model"
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! set atmospheric primitive variables !
prim2_a(:) = 0.0D0
eps2_a = 0.0D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

endsubroutine
