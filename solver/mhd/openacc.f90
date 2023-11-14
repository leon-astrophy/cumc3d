!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine handles data transfer between host and devices 
! AT THE BEGINNING of the simulations 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE POPULATE_DEVICE
USE MHD_MODULE
USE DEFINITION
USE RIEMANN_MODULE
IMPLICIT NONE

! Now populate all necessary, and reuseable arrays to the graphic cards !
!$ACC enter DATA COPYIN(boundary_flag, bfac_xin, bfac_yin, bfac_zin, bfac_xout, bfac_yout, bfac_zout, &
!$ACC x, y, z, xF, yF, zF, dx, dy, dz, xbar, dx_cb, dsine, dcose, vol, sine, sinf, wx, wy, wz, &
!$ACC prim, cons, cs, dpdrho, dpdeps, epsilon, sc, flux, l_rk, u_old, &
!$ACC eface, ecell, bcell, efield_x, efield_y, efield_z, &
!$ACC u_hll, p_hll, ustarL, ustarR, usstarL, usstarR, pstarL, pstarR, psstarL, psstarR, epsR, epsL, csL, csR, &
!$ACC primL, primR, fluxL, fluxR, uL, uR)

!$ACC UPDATE DEVICE(ggas, kgas)

! Populate arrays according to model !
CALL CUSTOM_POPULATE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine handles data transfer between host and devices 
! AT THE END of the simulations. We clear memories in the GPU
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CLEAR_DEVICE
USE MHD_MODULE
USE DEFINITION
USE RIEMANN_MODULE
IMPLICIT NONE
  
! Now we clear memory in the GPU device !
!$ACC exit DATA DELETE(boundary_flag, bfac_xin, bfac_yin, bfac_zin, bfac_xout, bfac_yout, bfac_zout, &
!$ACC x, y, z, xF, yF, zF, dx, dy, dz, xbar, dx_cb, dsine, dcose, vol, sine, sinf, wx, wy, wz, &
!$ACC prim, cons, cs, dpdrho, dpdeps, epsilon, sc, flux, l_rk, u_old, &
!$ACC eface, ecell, bcell, efield_x, efield_y, efield_z, &
!$ACC u_hll, p_hll, ustarL, ustarR, usstarL, usstarR, pstarL, pstarR, psstarL, psstarR, epsR, epsL, csL, csR, &
!$ACC primL, primR, fluxL, fluxR, uL, uR)

! Clear arrays according to model !
CALL CUSTOM_CLEAR

END SUBROUTINE