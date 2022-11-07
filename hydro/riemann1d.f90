!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Riemann problem tests for 1D hydrodynamics 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Riemann_1d
USE DEFINITION
IMPLICIT NONE

! Integer parameter !
integer :: j

! Set the polytropic index !
ggas2 = 1.4E0_DP

! Select according to models !
if(test_model == 1) then
	do j=1,nx_2
		if(j <= nx_2*3/10) then
			prim2(irho2,j,:,:) = 1.0E0_DP
			prim2(ivel2_x,j,:,:) = 0.75E0_DP
			prim2(itau2,j,:,:) = 1.0E0_DP
			epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2-1.0E0_DP)
		else
			prim2(irho2,j,:,:) = 0.125E0_DP
			prim2(ivel2_x,j,:,:) = 0.0E0_DP
			prim2(itau2,j,:,:) = 0.1E0_DP
			epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2-1.0E0_DP)
      		endif
	enddo
elseif(test_model == 2) then
	do j=1,nx_2
		if(j <= nx_2*5/10) then
			prim2(irho2,j,:,:) = 1.0E0_DP
			prim2(ivel2_x,j,:,:) = -2.0E0_DP
			prim2(itau2,j,:,:) = 0.4E0_DP
			epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2-1.0E0_DP)
		else
			prim2(irho2,j,:,:) = 1.0E0_DP
			prim2(ivel2_x,j,:,:) = 2.0E0_DP
			prim2(itau2,j,:,:) = 0.4E0_DP
			epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2-1.0E0_DP)
		endif
	enddo
elseif(test_model == 4) then
	do j=1,nx_2
		if(j <= nx_2*4/10) then
			prim2(irho2,j,:,:) = 5.99924E0_DP
			prim2(ivel2_x,j,:,:) = 19.5975E0_DP
			prim2(itau2,j,:,:) = 460.894E0_DP
			epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2-1.0E0_DP)
		else
			prim2(irho2,j,:,:) = 5.99242E0_DP
			prim2(ivel2_x,j,:,:) = -6.19633E0_DP
			prim2(itau2,j,:,:) = 46.0950E0_DP
			epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2-1.0E0_DP)
		endif
	enddo
elseif(test_model == 5) then
	do j=1,nx_2
		if(j <= nx_2*8/10) then
			prim2(irho2,j,:,:) = 1.0E0_DP
			prim2(ivel2_x,j,:,:) = -19.59745E0_DP
			prim2(itau2,j,:,:) = 1000.0E0_DP
			epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2-1.0E0_DP)
		else
			prim2(irho2,j,:,:) = 1.0E0_DP
			prim2(ivel2_x,j,:,:) = -19.59745E0_DP
			prim2(itau2,j,:,:) = 0.01E0_DP
			epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2-1.0E0_DP)
		endif
	enddo
elseif(test_model == 3) then
	do j=1,nx_2
		if(j <= nx_2*5/10) then
			prim2(irho2,j,:,:) = 1.0E0_DP
			prim2(ivel2_x,j,:,:) = 0.0E0_DP
			prim2(itau2,j,:,:) = 1000.0E0_DP
			epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2-1.0E0_DP)
		else
			prim2(irho2,j,:,:) = 1.0E0_DP
			prim2(ivel2_x,j,:,:) = 0.0E0_DP
			prim2(itau2,j,:,:) = 0.01E0_DP
			epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2-1.0E0_DP)
		endif
	enddo
elseif(test_model == 6) then
	do j=1,nx_2
		if(j <= nx_2*5/10) then
			prim2(irho2,j,:,:) = 1.4E0_DP
			prim2(ivel2_x,j,:,:) = 0.0E0_DP
			prim2(itau2,j,:,:) = 1.0E0_DP
			epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2-1.0E0_DP)
		else
			prim2(irho2,j,:,:) = 1.0E0_DP
			prim2(ivel2_x,j,:,:) = 0.0E0_DP
			prim2(itau2,j,:,:) = 1.0E0_DP
			epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2-1.0E0_DP)
		endif
	enddo
elseif(test_model == 7) then
	do j=1,nx_2
		if(j <= nx_2*5/10) then
			prim2(irho2,j,:,:) = 1.4E0_DP
			prim2(ivel2_x,j,:,:) = 0.1E0_DP
			prim2(itau2,j,:,:) = 1.0E0_DP
			epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2-1.0E0_DP)
		else
			prim2(irho2,j,:,:) = 1.0E0_DP
			prim2(ivel2_x,j,:,:) = 0.1E0_DP
			prim2(itau2,j,:,:) = 1.0E0_DP
			epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2-1.0E0_DP)
		endif
	enddo
elseif(test_model == 8) then
	do j=1,nx_2
		if(j <= nx_2*1/10) then
			prim2(irho2,j,:,:) = 1.0E0_DP
			prim2(ivel2_x,j,:,:) = 0.0E0_DP
			prim2(itau2,j,:,:) = 1000.0E0_DP
			epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2-1.0E0_DP)
		elseif(j >= nx_2*9/10) then
			prim2(irho2,j,:,:) = 1.0E0_DP
			prim2(ivel2_x,j,:,:) = 0.0E0_DP
			prim2(itau2,j,:,:) = 100.0E0_DP
			epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2-1.0E0_DP)
		else
			prim2(irho2,j,:,:) = 1.0E0_DP
			prim2(ivel2_x,j,:,:) = 0.0E0_DP
			prim2(itau2,j,:,:) = 0.01E0_DP
			epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2-1.0E0_DP)
		endif
	enddo
elseif(test_model == 9) then
	do j=1,nx_2
		if(j <= 4) then
			prim2(irho2,j,:,:) = 1.0E0_DP
			prim2(ivel2_x,j,:,:) = 0.0E0_DP
			prim2(itau2,j,:,:) = (3.0D0 * (ggas2 - 1.0D0) * 1.0D0)/ &
				(4.0D0 * pi * x2(4) ** 3)
			epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2-1.0E0_DP)
		else
			prim2(irho2,j,:,:) = 1.0E0_DP
			prim2(ivel2_x,j,:,:) = 0.0E0_DP
			prim2(itau2,j,:,:) = 1.0E-5_DP
			epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2-1.0E0_DP)
		endif
	enddo
elseif(test_model == 10) then
	do j=1,nx_2
		if(j < nx_2/8) then
			prim2(irho2,j,:,:) = 3.857143D0
			prim2(ivel2_x,j,:,:) = 2.629369D0
			prim2(itau2,j,:,:) = 31.0D0/3.0D0
			epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2-1.0E0_DP)
		else
			prim2(irho2,j,:,:) = 1.0E0_DP + 0.2D0*sin(8.0D0*(DBLE(j) - 0.5D0)*dx2)
			prim2(ivel2_x,j,:,:) = 0.0E0_DP
			prim2(itau2,j,:,:) = 1.0D0
			epsilon2(j,:,:) = prim2(itau2,j,:,:) / prim2(irho2,j,:,:) / (ggas2-1.0E0_DP)
		endif
	enddo
else
	stop "no such test model"
END IF

! set boundary conditions !
call BOUNDARYP_NM
call BOUNDARY1D_NM (epsilon2,even,part)

! set atmospheric primitive variables !
prim2_a(:) = 0.0D0
eps2_a = 0.0D0
temp2_a = 0.0D0

endsubroutine