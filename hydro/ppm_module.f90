!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the PPM (Piecewise Parabolic Method) Module that reconstruct   !
! interface values of primitive variables. We can choose either extremum !
! preserving limiter or the original limiter. We use shock flattening	 !
! and contact steepening algorithm. For details, please refer to the	 !
! original PPM paper and the method paper by FLASH code. Also see do	 !
! refer to the PPM reconstruction procedure by GR1D code		 !     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE PPM_MODULE
USE DEFINITION
IMPLICIT NONE

contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate cell average values to interface values using PPM 
! Assumed non-uniform gridding. Using the original PPM algorithm
! See Colella 1984. No steepening or flattening is performed.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PPM_1984 (i_in, grid_in, vm2, vm1, vc, vp1, vp2, vm_out, vp_out, type_in)
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER, INTENT(IN) :: i_in, type_in, grid_in

! The input into the subroutine, including conservative variable and input flux function !
REAL (DP), INTENT (IN) :: vm2, vm1, vc, vp1, vp2

! The output of the subroutine, the flux at cell boundary !
REAL (DP), INTENT (OUT) :: vp_out, vm_out

! Grid size !
REAL (DP) :: dm2, dm1, dc, dp1, dp2

! Temporal variables !
REAL (DP) :: vl, vr
REAL (DP) :: dmp1, dmc, dmm1
REAL (DP) :: deltap1, deltac, deltam1
REAL (DP) :: a1R, a2R, deltaXR, z1R, z2R
REAL (DP) :: a1L, a2L, deltaXL, z1L, z2L
REAL (DP) :: cp1, cp2, cp3
REAL (DP) :: cc1, cc2, cc3
REAL (DP) :: cm1, cm2, cm3
REAL (DP) :: condition

! Assign grid !
IF(type_in == 1) THEN
	IF(grid_in == 1) THEN
		dm2 = dr1(i_in - 2)
		dm1 = dr1(i_in - 1)
		dc = dr1(i_in)
		dp1 = dr1(i_in + 1)
		dp2 = dr1(i_in + 2)
	ELSEIF(grid_in == 2) THEN
		dm2 = dy1
		dm1 = dy1
		dc = dy1
		dp1 = dy1
		dp2 = dy1
	ELSEIF(grid_in == 3) THEN
		dm2 = dz1
		dm1 = dz1
		dc = dz1
		dp1 = dz1
		dp2 = dz1
	END IF
ELSEIF(type_in == 2) THEN
	IF(grid_in == 1) THEN
		dm2 = dr2(i_in - 2)
		dm1 = dr2(i_in - 1)
		dc = dr2(i_in)
		dp1 = dr2(i_in + 1)
		dp2 = dr2(i_in + 2)
	ELSEIF(grid_in == 2) THEN
		dm2 = dy2
		dm1 = dy2
		dc = dy2
		dp1 = dy2
		dp2 = dy2
	ELSEIF(grid_in == 3) THEN
		dm2 = dz2
		dm1 = dz2
		dc = dz2
		dp1 = dz2
		dp2 = dz2
	END IF
END IF

! Get coefficient !
a1R = dc/(dc + dp1)
a1L = dm1/(dm1 + dc)
a2R = 2.0D0*dp1*dc/(dp1 + dc)
a2L = 2.0D0*dc*dm1/(dc + dm1)
deltaXR = dm1 + dc + dp1 + dp2
deltaXL = dm2 + dm1 + dc + dp1
z1R = (dm1 + dc)/(2.0D0*dc + dp1)
z1L = (dm2 + dm1)/(2.0D0*dm1 + dc)
z2R = (dp2 + dp1)/(2.0D0*dp1 + dc)
z2L = (dp1 + dc)/(2.0D0*dc + dm1)

! For slope estimations !
cc1 = dc/(dm1 + dc + dp1)
cp1 = dp1/(dc + dp1 + dp2)
cm1 = dm1/(dm2 + dm1 + dc)
cc2 = (2.0D0*dm1 + dc)/(dp1 + dc)
cp2 = (2.0D0*dc + dp1)/(dp2 + dp1)
cm2 = (2.0D0*dm2 + dm1)/(dc + dm1)
cc3 = (dc + 2.0D0*dp1)/(dm1 + dc)
cp3 = (dp1 + 2.0D0*dp2)/(dc + dp1)
cm3 = (dm1 + 2.0D0*dc)/(dm2 + dm1)

! Get slopes !
deltap1 = cp1*(cp2*(vp2 - vp1) + cp3*(vp1 - vc))
deltac = cc1*(cc2*(vp1 - vc) + cc3*(vc - vm1))
deltam1 = cm1*(cm2*(vc - vm1) + cm3*(vm1 - vm2))

! Van Leer limiter !
condition = (vp1 - vc)*(vc - vm1)
IF(condition > 0.0D0) THEN
	dmc = min(abs(deltac), 2.0D0*abs(vc - vm1), 2.0D0*abs(vc - vp1))*sign(1.0D0,deltac)
ELSE
	dmc = 0.0D0
END IF

condition = (vp2 - vp1)*(vp1 - vc)
IF(condition > 0.0D0) THEN
	dmp1 = min(abs(deltap1), 2.0D0*abs(vp1 - vc), 2.0D0*abs(vp1 - vp2))*sign(1.0D0,deltap1)
ELSE
	dmp1 = 0.0D0
END IF

condition = (vc - vm1)*(vm1 - vm2)
IF(condition > 0.0D0) THEN
	dmm1 = min(abs(deltam1), 2.0D0*abs(vm1 - vm2), 2.0D0*abs(vm1 - vc))*sign(1.0D0,deltam1)
ELSE
	dmm1 = 0.0D0
END IF

! Get the interpolant !
vp_out = vc + a1R*(vp1 - vc) + (a2R*(z1R - z2R)*(vp1 - vc) - dc*z1R*dmp1 + dp1*z2R*dmc)/deltaXR
vm_out = vm1 + a1L*(vc - vm1) + (a2L*(z1L - z2L)*(vc - vm1) - dm1*z1L*dmc + dc*z2L*dmm1)/deltaXL

! backup !
vl = vm_out
vr = vp_out

! Check for extremum conditions !
IF((vr - vc)*(vc - vl) <= 0.0D0) THEN
	vp_out = vc
	vm_out = vc
ELSE
	condition = (vr - vl)*(vl - 3.0D0*vc + 2.0D0*vr)
	IF(condition < 0.0D0) THEN
		vm_out = 3.0d0*vc - 2.0D0*vr
	END IF
	condition = (vp_out - vm_out)*(3.0D0*vc - 2.0d0*vl - vr)
	IF(condition < 0.0D0) THEN
		vp_out = 3.0d0*vc - 2.0D0*vl
	END IF
END IF

END SUBROUTINE

END MODULE
