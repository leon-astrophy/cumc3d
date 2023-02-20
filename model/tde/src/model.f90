!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Builiding initial model 
! Newtonian torus with the Paczynsky-Wiita pseudo-Newtonean gravitational potential
! Solve the integral equations H + phi + h**2/(2-2q)/s**(2-2q) = c_const
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GET_MODEL
USE DEFINITION
IMPLICIT NONE
INCLUDE "param.h"

! Integer !
INTEGER :: i, j, k, l

! Real !
REAL*8 :: rad_sp
REAL*8 :: omega

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! override adiabatic index !
ggas2 = gamma

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! solve for the density profile !
DO l = 1, nz_2
  DO k = 1, ny_2
    DO j = 1, nx_2

      ! radius with respect to the star center
      rad_sp = (x_0*COS(y2(k)) - x2(j))**2 + (x_0*SIN(y2(k)))**2 + z2(l)**2
      rad_sp = SQRT(rad_sp)
      prim2(irho2,j,k,l) = rho_max*rad_sp**(-q_grad)
      omega = 1.0D0/(SQRT(x_0)*(x_0 - r_sh))
      prim2(ivel2_y,j,k,l) = omega*x_0*0.9D0

      ! Set density floor !
      IF(prim2(irho2,j,k,l) < rho_fac) THEN
        prim2(irho2,j,k,l) = rho_fac
      END IF

    END DO
  END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Find pressure and epsilon
DO l = 1, nz_2
  DO k = 1, ny_2
    DO j = 1, nx_2
      prim2(itau2,j,k,l) = p_max*(prim2(irho2,j,k,l)/rho_max)**(ggas2)
      epsilon2(j,k,l) = prim2(itau2,j,k,l)/prim2(irho2,j,k,l)/(ggas2 - 1.0d0)
    END DO
  END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! set atmospheric primitive variables !
prim2_a(:) = 0.0D0
prim2_a(irho2) = rho_fac
prim2_a(itau2) = p_max*(rho_fac/rho_max)**(ggas2)
eps2_a = prim2_a(itau2)/prim2_a(irho2)/(ggas2 - 1.0d0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! set boundary conditions !
call BOUNDARY1D_NM (epsilon2,part,even, even, even, even, even, even)
call BOUNDARYP_NM

! Set output profile interval !
total_time = 10.0d0
output_profiletime = 1.0d0

END SUBROUTINE
