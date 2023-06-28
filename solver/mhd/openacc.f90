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
!$ACC x2, y2, z2, xF2, yF2, zF2, dx2, dy2, dz2, x2cen, y2cen, x2bar, dx2_sq, dx2_cb, dsin2, dcos2, vol2, &
!$ACC prim2_a, prim2, cons2, cs2, dpdrho2, dpdeps2, epsilon2, sc2, flux_2, l2, u_old2, &
!$ACC lxm2, lxm1, lxc, lxp1, lxp2, lym2, lym1, lyc, lyp1, lyp2, lzm2, lzm1, lzc, lzp1, lzp2, &
!$ACC rxm2, rxm1, rxc, rxp1, rxp2, rym2, rym1, ryc, ryp1, ryp2, rzm2, rzm1, rzc, rzp1, rzp2, &
!$ACC hpx, hpy, hpz, hmx, hmy, hmz, &
!$ACC eface, ecell_x, ecell_y, ecell_z, bcell, efield_x, efield_y, efield_z, mflux_x, &
!$ACC u_hll, p_hll, ustarL, ustarR, usstarL, usstarR, pstarL, pstarR, psstarL, psstarR, eps2R, eps2L, cs2L, cs2R, &
!$ACC primL2, primR2, fluxL2, fluxR2, uL2, uR2)

!$ACC UPDATE DEVICE(ggas2, kgas2)

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
!$ACC x2, y2, z2, xF2, yF2, zF2, dx2, dy2, dz2, x2cen, y2cen, x2bar, dx2_sq, dx2_cb, dsin2, dcos2, vol2, &
!$ACC prim2_a, prim2, cons2, cs2, dpdrho2, dpdeps2, epsilon2, sc2, flux_2, l2, u_old2, &
!$ACC lxm2, lxm1, lxc, lxp1, lxp2, lym2, lym1, lyc, lyp1, lyp2, lzm2, lzm1, lzc, lzp1, lzp2, &
!$ACC rxm2, rxm1, rxc, rxp1, rxp2, rym2, rym1, ryc, ryp1, ryp2, rzm2, rzm1, rzc, rzp1, rzp2, &
!$ACC hpx, hpy, hpz, hmx, hmy, hmz, &
!$ACC eface, ecell_x, ecell_y, ecell_z, bcell, efield_x, efield_y, efield_z, mflux_x, &
!$ACC u_hll, p_hll, ustarL, ustarR, usstarL, usstarR, pstarL, pstarR, psstarL, psstarR, eps2R, eps2L, cs2L, cs2R, &
!$ACC primL2, primR2, fluxL2, fluxR2, uL2, uR2)

! Clear arrays according to model !
CALL CUSTOM_CLEAR

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine handles update devices memory from the host
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DATA_HOST2DEVICE(array_in)
USE MHD_MODULE
USE DEFINITION
USE RIEMANN_MODULE
IMPLICIT NONE

! Define target array !
REAL*8, INTENT (IN), DIMENSION (-2:nx_2+3,-2:ny_2+3,-2:nz_2+3) :: array_in

! Now do the data transfer !
!$ACC UPDATE DEVICE(array_in)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine handles update host memory from the device
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DATA_DEVICE2HOST(array_in)
USE MHD_MODULE
USE DEFINITION
USE RIEMANN_MODULE
IMPLICIT NONE
  
! Define target array !
REAL*8, INTENT (IN), DIMENSION (-2:nx_2+3,-2:ny_2+3,-2:nz_2+3) :: array_in

! Now do the data transfer !
!$ACC UPDATE HOST(array_in)

END SUBROUTINE